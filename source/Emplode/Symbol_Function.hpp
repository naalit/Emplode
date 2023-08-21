/**
 *  @note This file is part of Emplode, currently within https://github.com/mercere99/MABE2
 *  @copyright Copyright (C) Michigan State University, MIT Software license; see doc/LICENSE.md
 *  @date 2019-2021.
 *
 *  @file  Symbol_Function.hpp
 *  @brief Manages individual functions for config.
 *  @note Status: BETA
 */

#ifndef EMPLODE_SYMBOL_FUNCTION_HPP
#define EMPLODE_SYMBOL_FUNCTION_HPP

#include <string>

#include "AST.hpp"
#include "emp/base/Ptr.hpp"
#include "emp/base/notify.hpp"
#include "emp/base/vector.hpp"
#include "emp/datastructs/tuple_utils.hpp"
#include "emp/meta/ValPack.hpp"

#include "Symbol.hpp"
#include "SymbolTableBase.hpp"

namespace emplode {

  class Symbol_Function : public Symbol {
  protected:
    using symbol_ptr_t = emp::Ptr<Symbol>;
    using fun_t = symbol_ptr_t( const emp::vector<symbol_ptr_t> & );
    using std_fun_t = std::function< fun_t >;

  private:
    using this_t = Symbol_Function;

    struct FunInfo {
      std_fun_t fun;            // Unified-form function.
      int num_params;        // How many arguments does this function take?
    };

    emp::vector<FunInfo> overloads;  // Set of overload options for this function.
    emp::TypeID return_type;         // All overloads must share a return type.

    // size_t arg_count;

  public:
    Symbol_Function(const std::string & _name,
                    std_fun_t fun,
                    const std::string & _desc,
                    emp::Ptr<Symbol_Scope> _scope,
                    int num_params,
                    emp::TypeID _ret_type)
      : Symbol(_name, _desc, _scope), return_type(_ret_type)
    {
      overloads.push_back( FunInfo{fun, num_params} );
    }

    Symbol_Function(const Symbol_Function &) = default;
    emp::Ptr<Symbol> Clone() const override { return emp::NewPtr<this_t>(*this); }

    std::string GetTypename() const override { return "[Symbol_Function]"; }

    bool IsFunction() const override { return true; }
    bool HasNumericReturn() const override { return return_type.IsArithmetic(); }
    bool HasStringReturn() const override { return return_type.IsType<std::string>(); }

    std::string AsString() const override { return "[[__FUNCTION__]]"; }

    /// Set this symbol to be a correctly-typed scope pointer.
    emp::Ptr<Symbol_Function> AsFunctionPtr() override { return this; }
    emp::Ptr<const Symbol_Function> AsFunctionPtr() const override { return this; }

    bool CopyValue(const Symbol & in) override {
      if (in.IsFunction() == false) {
          std::cerr << "Trying to assign `" << in.GetName() << "' to '" << GetName()
                    << "', but " << in.GetName() << " is not a Function." << std::endl;
        return false;   // Mis-matched types; failed to copy.
      }

      const Symbol_Function & in_fun = in.AsFunction();
      overloads = in_fun.overloads;

      return true;
    }


    symbol_ptr_t Call( const emp::vector<symbol_ptr_t> & args ) override {
      emp_assert(overloads.size() > 0);

      // Find the correct overloads...
      for (auto & x : overloads) {
        if (x.num_params == -1 || x.num_params == (int) args.size()) return x.fun(args);
      }

      std::string msg =
        emp::to_string("No overload for function '", name, "' that takes ", args.size(),
                       " arguments.\n...", overloads.size(), " options are:");
      for (const auto & x : overloads) {
        msg += emp::to_string (' ', x.num_params);
      }
      emp::notify::Exception("mabe::Symbol_Function::NO_OVERLOAD", msg);
      return nullptr;
    }
  };

  class Symbol_UserFunction : public Symbol_Function {
  private:
    emp::Ptr<Symbol_Scope> own_scope;
    emp::Ptr<ASTNode_Block> body;

    static std_fun_t make_fun(emp::Ptr<ASTNode_Block> body, emp::vector<Var> &params) {
      return [body, params](const emp::vector<emp::Ptr<Symbol>> & args) {
          if (args.size() != params.size()) {
            std::cerr << "Expected " << params.size() << " arguments but got " << args.size() << std::endl;
            exit(1);
          }
          for (int i = 0; i < args.size(); i++) {
            auto param = params[i];
            param.SetValue(args[i]->ShallowClone());
          }

          auto result = body->Process();
          if (result && result->IsReturn()) {
            auto ret = result.DynamicCast<Symbol_Special>()->ReturnValue();
            result.Delete();
            return ret;
          } else {
            emp::Ptr<Symbol> ret = nullptr;
            return ret;
          }
        };
    }

  public:
    Symbol_UserFunction(const std::string & _name,
                    emp::vector<Var> params,
                    emp::Ptr<ASTNode_Block> body,
                    const std::string & _desc,
                    emp::Ptr<Symbol_Scope> _scope,
                    emp::TypeID _ret_type,
                    emp::Ptr<Symbol_Scope> own_scope)
      : own_scope(own_scope), body(body),
      Symbol_Function(_name, make_fun(body, params), _desc, _scope, params.size(), _ret_type) {}

    ~Symbol_UserFunction() {
      own_scope.Delete();
      body.Delete();
    }
  };

}

#endif
