/**
 *  @note This file is part of Emplode, currently within https://github.com/mercere99/MABE2
 *  @copyright Copyright (C) Michigan State University, MIT Software license; see doc/LICENSE.md
 *  @date 2019-2021.
 *
 *  @file  Symbol_Scope.hpp
 *  @brief Manages a full scope with many internal symbols (including sub-scopes).
 *  @note Status: BETA
 * 
 *  DEVELOPER NOTES:
 *  - Need to fix Add() function to give a user-level error, rather than an assert on duplication.
 */

#ifndef EMPLODE_SYMBOL_SCOPE_HPP
#define EMPLODE_SYMBOL_SCOPE_HPP

#include "AST.hpp"
#include "emp/base/map.hpp"
#include "emp/datastructs/map_utils.hpp"
#include <optional>

#include "Symbol.hpp"
#include "Symbol_Function.hpp"
#include "Symbol_Linked.hpp"
#include "TypeInfo.hpp"

namespace emplode {

  class EmplodeType;
  class Symbol_Object;

  // Set of multiple config entries.
  class Symbol_Scope : public Symbol {
  protected:
    using symbol_ptr_t = emp::Ptr<Symbol>;
    using const_symbol_ptr_t = emp::Ptr<const Symbol>;
    emp::map< std::string, Var > symbol_map;   ///< Map of names to entries.

    template <typename T, typename... ARGS>
    Var Add(const std::string & name, ARGS &&... args) {
      auto new_ptr = emp::NewPtr<T>(name, std::forward<ARGS>(args)...);
      emp_always_assert(!emp::Has(symbol_map, name), "Do not redeclare functions or variables!",
                 name);
      auto entry = symbol_map.insert({name, Var(new_ptr)});
      return entry.first->second;
    }

    template <typename T, typename... ARGS>
    Var AddBuiltin(const std::string & name, ARGS &&... args) {
      Var result = Add<T>(name, std::forward<ARGS>(args)...);
      result.GetValue()->SetBuiltin();
      return result;
    }

  public:
    Symbol_Scope(const std::string & _name, const std::string & _desc, emp::Ptr<Symbol_Scope> _scope)
      : Symbol(_name, _desc, _scope) { }

    std::string GetTypename() const override { return "Scope"; }

    bool IsScope() const override { return true; }
    bool IsLocal() const override { return true; }  // @CAO, for now assuming all scopes are local!

    std::string AsString() const override { return "[[__SCOPE__]]"; }

    /// Set this symbol to be a correctly-typed scope pointer.
    emp::Ptr<Symbol_Scope> AsScopePtr() override { return this; }
    emp::Ptr<const Symbol_Scope> AsScopePtr() const override { return this; }

    bool CopyValue(const Symbol & in) override {
      if (in.IsScope() == false) {
          std::cerr << "Trying to assign `" << in.GetName() << "' to '" << GetName()
                    << "', but " << in.GetName() << " is not a Scope." << std::endl;
        return false;   // Mis-matched types; failed to copy.
      }

      const Symbol_Scope & in_scope = in.AsScope();

      // Assignment to an existing Struct cannot create new variables; all must already exist.
      // Do not delete other existing entries.
      for (const auto & [name, var] : in_scope.symbol_map) {
        // If entry does not exist fail the copy.
        if (!emp::Has(symbol_map, name)) {
          std::cerr << "Trying to assign `" << in.GetName() << "' to '" << GetName()
                    << "', but " << GetName() << "." << name << " does not exist." << std::endl;
          return false;
        }

        if (var.GetValue()->IsFunction()) continue; // Don't copy functions.

        bool success = GetSymbol(name)->GetValue()->CopyValue(*var.GetValue());
        if (!success) {
          std::cerr << "Trying to assign `" << in.GetName() << "' to '" << GetName()
                    << "', but failed on `" << GetName() << "." << name << "`." << std::endl;
          return false; // Stop immediately on failure.
        }
      }

      // If we made it this far, it must have worked!
      return true;
    }

    void CopyFields(const Symbol_Scope & in) {
      for (const auto & [name, var] : in.symbol_map) {
        if (symbol_map.contains(name)) {
          symbol_map.at(name).SetValue(var.GetValue()->ShallowClone());
        } else {
          if (var.GetValue()->GetDesc() == "Local variable created by assignment") {
            std::cerr << "Assignment to nonexistent member '" << name << "' of object or struct when initializing" << std::endl;
            exit(1);
          }
          symbol_map.insert({name, Var(var.GetValue()->ShallowClone())});
        }
      }
    }

    /// Get a symbol out of this scope; 
    std::optional<Var> GetSymbol(std::string name) {
      auto result = symbol_map.find(name);
      if (result == symbol_map.end()) {
        return {};
      } else {
        return result->second;
      }
    }

    /// Lookup a variable, scanning outer scopes if needed
    std::optional<Var> LookupSymbol(const std::string & name, bool scan_scopes=true) const {
      // See if this next symbol is in the var list.
      auto it = symbol_map.find(name);

      // If this name is unknown, check with the parent scope!
      if (it == symbol_map.end()) {
        if (scope.IsNull() || !scan_scopes) return {};  // No parent?  Just fail...
        return scope->LookupSymbol(name);
      }

      // Otherwise we found it!
      return it->second;
    }

    /// Add a configuration symbol that is linked to a variable - the incoming variable sets
    /// the default and is automatically updated when configs are loaded.
    template <typename VAR_T>
    Var LinkVar(const std::string & name,
                                        VAR_T & var,
                                        const std::string & desc,
                                        bool is_builtin = false) {
      if (is_builtin) return AddBuiltin<Symbol_Linked<VAR_T>>(name, var, desc, this);
      return Add<Symbol_Linked<VAR_T>>(name, var, desc, this);
    }

    /// Add a configuration symbol that interacts through a pair of functions - the functions are
    /// automatically called any time the symbol value is accessed (get_fun) or changed (set_fun)
    template <typename VAR_T>
    Var LinkFuns(const std::string & name,
                                            std::function<VAR_T()> get_fun,
                                            std::function<void(const VAR_T &)> set_fun,
                                            const std::string & desc,
                                            bool is_builtin = false) {
      if (is_builtin) {
        return AddBuiltin<Symbol_LinkedFunctions<VAR_T>>(name, get_fun, set_fun, desc, this);
      }
      return Add<Symbol_LinkedFunctions<VAR_T>>(name, get_fun, set_fun, desc, this);
    }

    /// Add an internal variable of type String.
    Var AddLocalVar(const std::string & name, const std::string & desc) {
      return Add<Symbol_Var>(name, 0.0, desc, this);
    }

    /// Add an internal scope inside of this one.
    Var AddScope(const std::string & name, const std::string & desc) {
      return Add<Symbol_Scope>(name, desc, this);
    }

    /// Add an internal scope inside of this one (defined in Symbol_Object.hpp)
    Var AddObject(
      const std::string & name,
      const std::string & desc,
      emp::Ptr<EmplodeType> obj_ptr,
      TypeInfo & type_info,
      bool obj_owned=false
    );

    template <typename FUN_T>
    int constexpr CountParams() {
      using info_t = emp::FunInfo<FUN_T>;

      // If we have a single argument and it's a vector of symbols, assume variable args.
      if constexpr (info_t::num_args == 1) {
        using arg_t = typename info_t::template arg_t<0>;
        if constexpr (std::is_same_v<arg_t, const emp::vector<symbol_ptr_t> &>) {
          return -1;
        }             
      }

      return info_t::num_args;
    }

    /// Add a new user-defined function.
    template <typename FUN_T>
    Var AddFunction(const std::string & name,  FUN_T fun,
                                  const std::string & desc,  emp::TypeID return_type) {
      return Add<Symbol_Function>(name, fun, desc, this, CountParams<FUN_T>(), return_type);
    }

    /// Add a new user-defined function, providing number of params separately.
    Var AddUserFunction(const std::string & name, const std::string & desc,
                                  emp::TypeID return_type, emp::vector<Var> params,
                                  emp::Ptr<ASTNode_Block> body, emp::Ptr<Symbol_Scope> scope) {
      return Add<Symbol_UserFunction>(name, params, body, desc, this, return_type, scope);
    }

    /// Add a new function that is a standard part of the scripting language.
    template <typename FUN_T>
    Var AddBuiltinFunction(const std::string & name,  FUN_T fun,
                                         const std::string & desc,  emp::TypeID return_type) {
      return AddBuiltin<Symbol_Function>(name, fun, desc, this, CountParams<FUN_T>(), return_type);
    }

    /// Write out all of the parameters contained in this scope to the provided stream.
    const Symbol & WriteContents(std::ostream & os=std::cout, const std::string & prefix="",
                                      size_t comment_offset=32) const {

      // Loop through all of the entires in this scope and Write them.
      for (auto [name, var] : symbol_map) {
        if (var.GetValue()->IsBuiltin()) continue; // Skip writing built-in entries.
        var.GetValue()->Write(os, prefix, comment_offset);
      }

      return *this;
    }

    /// Write out this scope AND it's contents to the provided stream.
    const Symbol & Write(std::ostream & os=std::cout, const std::string & prefix="",
                         size_t comment_offset=32) const override
    {
      // If this is a built-in scope, don't print it.
      if (IsBuiltin()) return *this;

      // Declare this scope, starting with the type if originally declared locally.
      std::string cur_line = prefix;
      if (IsLocal()) cur_line += emp::to_string(GetTypename(), " ");
      cur_line += name;

      bool has_body = emp::AnyOf(symbol_map, [](Var var){ return !var.GetValue()->IsBuiltin(); });

      // Only open this scope if there are contents.
      cur_line += has_body ? " { " : ";";
      os << cur_line;

      // Indent the comment for the description (if there is one)
      WriteDesc(os, comment_offset, cur_line.size());

      // If we have internal entries, write them out.
      if (has_body) {
        WriteContents(os, prefix+"  ", comment_offset);
        os << prefix << "}\n";      // Close the scope.
      }

      return *this;
    }

    /// Make a copy of this scope and all of the entries inside it.
    symbol_ptr_t Clone() const override {
      emp::Ptr<Symbol_Scope> result = emp::NewPtr<Symbol_Scope>(name, desc, scope);
      for (auto [name, var] : symbol_map) {
        result->symbol_map.insert({name, Var(var.GetValue()->Clone())});
      }
      return result;
    }

    symbol_ptr_t ShallowClone() const override { return emp::NewPtr<Symbol_Scope>(*this); }
  };

  class Symbol_List : public Symbol {
  private:
    std::shared_ptr<emp::vector<symbol_ptr_t>> values;
    emp::Ptr<Symbol_Scope> member_funs;

    template<size_t arity, typename T>
    void AddMemberFun(std::string name, std::string desc, emp::TypeID ret_type, T fun) {
      auto wrapped = [name, fun](const std::vector<symbol_ptr_t> &input) {
        if (input.size() != arity) {
          std::cerr << "Error: list method '" << name << "' takes " << arity << " arguments, but was given " << input.size() << std::endl;
          exit(1);
        }
        if constexpr (arity == 0) return fun();
        if constexpr (arity == 1) return fun(input[0]);
        if constexpr (arity == 2) return fun(input[0], input[1]);
      };
      member_funs->AddBuiltinFunction(name, wrapped, desc, ret_type);
    }

  public:
    Symbol_List() : Symbol("__List", "List", nullptr) {
      values = std::make_shared<emp::vector<symbol_ptr_t>>();
      member_funs.New("List", "List scope", nullptr);
      AddMemberFun<1>("push", "Add a value to the end of the list", emp::GetTypeID<void>(), [this](symbol_ptr_t value) {
        if (!value->IsTemporary()) {
          value = value->ShallowClone();
        }
        value->SetTemporary(false);
        values->push_back(value);
        return nullptr;
      });
      AddMemberFun<0>("pop", "Removes and returns the value at the end of the list", emp::GetTypeID<symbol_ptr_t>(), [this]() {
        symbol_ptr_t value = values->back();
        value->SetTemporary();
        values->pop_back();
        return value;
      });
    }

    ~Symbol_List() {
      if (values.use_count() == 1) {
        for (auto i : *values) {
          i.Delete();
        }
      }
      member_funs.Delete();
    }

    std::string GetTypename() const override { return "List"; }

    emp::Ptr<Symbol_Scope> AsScopePtr() override {
      return member_funs;
    }

    symbol_ptr_t Clone() const override {
      auto list = emp::NewPtr<Symbol_List>();
      // TODO Should this call Clone() or ShallowClone() on the inner values?
      for (auto i : *values) {
        list->Push(i->Clone());
      }
      return list;
    }

    symbol_ptr_t ShallowClone() const override {
      auto list = emp::NewPtr<Symbol_List>();
      list->values = values;
      return list;
    }

    void Push(symbol_ptr_t value) {
      values->push_back(value);
    }

    LValue Get(size_t idx) {
      if (idx >= values->size()) {
        std::cerr << "index " << idx << " out of bounds for list of length " << values->size() << std::endl;
        exit(1);
      }

      return LValue(&(*values)[idx]);
    }

    void Print(std::ostream &os) const override {
      os << '[';
      bool first = true;
      for (auto i : *values) {
        if (!first)
          os << ", ";
        i->Print(os);
        first = false;
      }
      os << ']';
    }
  };

  // These has to be here (or in another downstream file) because of include cycle issues
  emp::Ptr<Symbol> ASTNode_ListInit::Process() {
    #ifndef NDEBUG
    emp::notify::Verbose(
      "Emplode::AST",
      "AST: Processing list initializer"
    );
    #endif
    auto list = emp::NewPtr<Symbol_List>();
    list->SetTemporary();
    for (auto i : children) {
      symbol_ptr_t value = i->Process();
      if (!value) {
        std::cerr << "Expression in list initializer does not produce a value" << std::endl;
        exit(1);
      }
      auto clone = value->ShallowClone();
      clone->SetTemporary(false);
      list->Push(clone);
      if (value->IsTemporary()) {
        value.Delete();
      }
    }
    return list;
  }

  emp::Ptr<Symbol> ASTNode_StructInit::Process() {
    #ifndef NDEBUG
    emp::notify::Verbose(
      "Emplode::AST",
      "AST: Processing struct initializer"
    );
    #endif
    auto scope = emp::NewPtr<Symbol_Scope>(name, "Local struct", nullptr);
    scope->SetTemporary();
    return scope;
  }

  emp::Ptr<Symbol> ASTNode_CopyFields::Process() {
    #ifndef NDEBUG
    emp::notify::Verbose(
      "Emplode::AST",
      "AST: Processing copy fields"
    );
    #endif
    emp_assert(children.size() == 2);
    symbol_ptr_t lhs = children[0]->Process();
    symbol_ptr_t rhs = children[1]->Process();
    auto lhs_scope = lhs->AsScopePtr();
    auto rhs_scope = rhs->AsScopePtr();
    if (!lhs_scope || !rhs_scope) {
      std::cerr << "ASTNode_CopyFields needs two scope operands" << std::endl;
      exit(1);
    }
    lhs_scope->CopyFields(*rhs_scope);
    // Rhs might be temporary and should be deleted, but lhs is the target and should never be temporary
    if (rhs->IsTemporary())
      rhs.Delete();
    return nullptr;
  }

  std::optional<LValue> ASTNode_Member::AsLValue() {
    emp_assert(children.size() == 1);

    auto scope = children[0]->Process()->AsScopePtr();
    if (!scope) {
      std::cerr << "tried to access member of non-scope:" << std::endl;
      PrintAST(std::cerr);
      exit(1);
    }
    // In C++23 this would use and_then()
    auto var = scope->GetSymbol(name);
    if (var.has_value())
      return var->AsLValue();
    else
      return {};
  }

  std::optional<LValue> ASTNode_Subscript::AsLValue() {
    #ifndef NDEBUG
    emp::notify::Verbose(
      "Emplode::AST",
      "AST: Processing subscript"
    );
    #endif

    emp_assert(children.size() == 2);

    auto list = children[0]->Process().DynamicCast<Symbol_List>();
    if (!list) {
      std::cerr << "tried to use subscript on non-list:" << std::endl;
      PrintAST(std::cerr);
      exit(1);
    }

    auto idx = children[1]->Process()->AsDouble();
    if (size_t(idx) != idx) {
      std::cerr << "tried to use subscript with negative or non-integer index " << idx << std::endl;
      PrintAST(std::cerr);
      exit(1);
    }

    return list->Get(idx);
  }
}
#endif
