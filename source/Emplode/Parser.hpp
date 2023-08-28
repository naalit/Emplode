/**
 *  @note This file is part of Emplode, currently within https://github.com/mercere99/MABE2
 *  @copyright Copyright (C) Michigan State University, MIT Software license; see doc/LICENSE.md
 *  @date 2021-2022.
 *
 *  @file  Parser.hpp
 *  @brief Manages parsing of Emplode language input streams.
 *  @note Status: BETA
 * 
 * 
 *  ParseState is a class that tracks:
 *  - The position in the incoming token stream that we are parsing.
 *  - The symbol table that is being constructed
 *  - The current scope we are in
 *  - The lexer that was used to build the token stream.
 * 
 *  The Parser itself uses ParseState to convert the incoming token stream into
 *  an abstract syntax tree.
 * 
 *  All Parse* functions take in the current parse state and advance it appropriately
 *  (including updating the symbol table and position in the token stream).  
 *  Unless otherwise specified, parse functions generate the abstract syntax tree nodes
 *  that represent the material parsed.  
 * 
 *  Parsing member functions include:
 * 
 *    ASTPtr ParseVar(ParseState & state, bool create_ok=false, bool scan_scopes=true);
 *    ASTPtr ParseValue(ParseState & state);
 *    ASTPtr ParseExpression(ParseState & state, bool decl_ok=false, size_t prec_limit=1000);
 *    Symbol & ParseDeclaration(ParseState & state);
 *    ASTPtr ParseEvent(ParseState & state);
 *    ASTPtr ParseKeywordStatement(ParseState & state);   // IF, WHILE, etc
 *    ASTPtr ParseStatement(ParseState & state);  // variable declaration, expression, or event.
 *    ASTPtr ParseStatementList(ParseState & state); // Go to end of scope of file
 * 
 *  Other helper functions include:
 *    /// Calculate the result of the provided operation on two computed entries.
 *    ASTPtr ProcessOperation(const emp::Token & op_token, ASTPtr tree1, ASTPtr tree2);
 */

#ifndef EMPLODE_PARSER_HPP
#define EMPLODE_PARSER_HPP

#include <string>
#include <utility>

#include "emp/base/Ptr.hpp"
#include "emp/base/vector.hpp"
#include "emp/tools/string_utils.hpp"

#include "AST.hpp"
#include "Lexer.hpp"
#include "Symbol_Scope.hpp"
#include "SymbolTable.hpp"

namespace emplode {

  class ParseState {
  private:
    emp::TokenStream::Iterator pos;
    emp::Ptr<SymbolTable> symbol_table;
    emp::vector< emp::Ptr<Symbol_Scope> > scope_stack;
    emp::Ptr<Lexer> lexer;

  public:
    ParseState(emp::TokenStream::Iterator _pos, SymbolTable & _table,
               Symbol_Scope & _scope, Lexer & _lexer)
      : pos(_pos), symbol_table(&_table), lexer(&_lexer) { scope_stack.push_back(&_scope); }
    ParseState(const ParseState &) = default;
    ~ParseState() { }

    ParseState & operator=(const ParseState &) = default;

    bool operator==(const ParseState & in) const { return pos == in.pos; }
    bool operator!=(const ParseState & in) const { return pos != in.pos; }
    bool operator< (const ParseState & in) const { return pos <  in.pos; }
    bool operator<=(const ParseState & in) const { return pos <= in.pos; }
    bool operator> (const ParseState & in) const { return pos >  in.pos; }
    bool operator>=(const ParseState & in) const { return pos >= in.pos; }

    ParseState & operator++() { ++pos; return *this; }
    ParseState operator++(int) { ParseState old(*this); ++pos; return old; }
    ParseState & operator--() { --pos; return *this; }
    ParseState operator--(int) { ParseState old(*this); --pos; return old; }

    bool IsValid() const { return pos.IsValid(); }
    bool AtEnd() const { return pos.AtEnd(); }

    size_t GetIndex() const { return pos.GetIndex(); }  ///< Return index in token stream.
    int GetLine() const { return (int) pos->line_id; }
    size_t GetTokenSize() const { return pos.IsValid() ? pos->lexeme.size() : 0; }
    SymbolTable & GetSymbolTable() { return *symbol_table; }
    Symbol_Scope & GetScope() {
      emp_assert(scope_stack.size() && scope_stack.back() != nullptr);
      return *scope_stack.back();
    }
    const Symbol_Scope & GetScope() const {
      emp_assert(scope_stack.size() && scope_stack.back() != nullptr);
      return *scope_stack.back();
    }
    const std::string & GetScopeName() const { return GetScope().GetName(); }

    std::string AsString() {
      return emp::to_string("[pos=", pos.GetIndex(),
                            ",lex='", AsLexeme(),
                            "',scope='", GetScope().GetName(),
                            "']");
    }

    bool IsKeyword() const { return pos && lexer->IsKeyword(*pos); }    
    bool IsID() const { return pos && lexer->IsID(*pos); }
    bool IsNumber() const { return pos && lexer->IsNumber(*pos); }
    bool IsString() const { return pos && lexer->IsString(*pos); }

    bool IsSignal() const { return symbol_table->HasSignal(AsLexeme()); }
    bool IsType() const { return symbol_table->HasType(AsLexeme()); }

    /// Convert the current state to a character; use \0 if cur token is not a symbol.
    char AsChar() const { return (pos && lexer->IsSymbol(*pos)) ? pos->lexeme[0] : 0; }

    /// Return the token associate with the current state.
    emp::Token AsToken() const { return *pos; }

    /// Return the lexeme associate with the current state.
    const std::string & AsLexeme() const { return pos ? pos->lexeme : emp::empty_string(); }

    /// Return the lexeme associate with the current state AND advance the token stream.
    const std::string & UseLexeme() {
      const std::string & out = AsLexeme();
      pos++;
      return out;
    }

    /// Return whether the current token is the specified lexeme; if so also advance token stream.
    bool UseIfLexeme(const std::string & test_str) {
      if (AsLexeme() != test_str) return false;
      ++pos;
      return true;
    }

    /// Return whether the current token is the specified char; if so also advance token stream.
    bool UseIfChar(char test_char) {
      if (AsChar() != test_char) return false;
      ++pos;
      return true;
    }

    /// Report an error in parsing this file and exit.
    template <typename... Ts>
    void Error(Ts &&... args) const {
      std::string line_info = pos.AtEnd() ? "end of input" : emp::to_string("line ", pos->line_id);
      
      emp::notify::Error("(", line_info, " in '", pos.GetTokenStream().GetName(), "'): ",
                         emp::to_string(std::forward<Ts>(args)...), "\nAborting.");
      exit(1);
    }

    template <typename... Ts>
    void Require(bool test, Ts &&... args) const {
      if (!test) Error(std::forward<Ts>(args)...);
    }
    template <typename... Ts>
    void RequireID(Ts &&... args) const {
      if (!IsID()) Error(std::forward<Ts>(args)...);
    }
    template <typename... Ts>
    void RequireNumber(Ts &&... args) const {
      if (!IsNumber()) Error(std::forward<Ts>(args)...);
    }
    template <typename... Ts>
    void RequireString(Ts &&... args) const {
      if (!IsString()) Error(std::forward<Ts>(args)...);
    }
    template <typename... Ts>
    void RequireChar(char req_char, Ts &&... args) const {
      if (AsChar() != req_char) Error(std::forward<Ts>(args)...);
    }
    template <typename... Ts>
    void RequireLexeme(const std::string & lex, Ts &&... args) const {
      if (AsLexeme() != lex) Error(std::forward<Ts>(args)...);
    }

    template <typename... Ts>
    void UseRequiredChar(char req_char, Ts &&... args) {
      if (AsChar() != req_char) Error(std::forward<Ts>(args)...);
      ++pos;      
    }

    void PushScope(Symbol_Scope & _scope) { scope_stack.push_back(&_scope); }
    void PopScope() { scope_stack.pop_back(); }

    Var LookupSymbol(const std::string & var_name, bool scan_scopes) {
      std::optional<Var> out_symbol = GetScope().LookupSymbol(var_name, scan_scopes);
      // If we can't find this identifier, throw an error.
      if (!out_symbol.has_value()) {
        Error("'", var_name, "' does not exist as a parameter, variable, or type.",
              "  Current scope is '", GetScope().GetName(), "'");
      }
      return out_symbol.value();
    }

    Var AddLocalVar(const std::string & name, const std::string & desc) {
      return GetScope().AddLocalVar(name, desc);
    }
    Var AddScope(const std::string & name, const std::string & desc) {
      return GetScope().AddScope(name, desc);
    }
    Var AddObject(const std::string & type_name, const std::string & var_name) {
      return symbol_table->MakeObjSymbol(type_name, var_name, GetScope());
    }

    /// Add a user-defined function
    Var AddFunction(const std::string & name, const std::string & desc,
        const std::string & ret_type, emp::vector<Var> params, emp::Ptr<ASTNode_Block> body,
        emp::Ptr<Symbol_Scope> scope) {
      emp::TypeID ret_type_id = symbol_table->GetType(ret_type).GetTypeID();
      return GetScope().AddUserFunction(name, desc, ret_type_id, params, body, scope);
    }

    /// Add an instance of an event with an action that should be triggered.
    template <typename... Ts>
    void AddAction(Ts &&... args) { symbol_table->AddAction(std::forward<Ts>(args)...); }
  };



  class Parser {
  private:
    std::unordered_map<std::string, size_t> precedence_map;  ///< Precedence levels for symbols.

    /// Print only when debugging.
    /// To activate debugging data, do: emp::notify::SetVerbose("emplode::Parser");
    template <typename... Ts>
    void Debug(Ts... args) const {
      emp::notify::Verbose("Emplode::Parser", emp::to_string(std::forward<Ts>(args)...));
    }

  public:
    Parser() {
      // Setup operator precedence.
      size_t cur_prec = 0;
      precedence_map["("] = precedence_map["["] = cur_prec++;
      precedence_map["**"] = cur_prec++;
      precedence_map["*"] = precedence_map["/"] = precedence_map["%"] = cur_prec++;
      precedence_map["+"] = precedence_map["-"] = cur_prec++;
      precedence_map["<"] = precedence_map["<="] = precedence_map[">"] = precedence_map[">="] = cur_prec++;
      precedence_map["=="] = precedence_map["!="] = cur_prec++;
      precedence_map["&&"] = cur_prec++;
      precedence_map["||"] = cur_prec++;
      precedence_map["="] = cur_prec++;
    }
    ~Parser() {}

    /// Load a variable name from the provided scope.
    /// @param state the current state of parsing (token stream, symbol table, etc.)
    /// @param create_ok indicates if we should create any variables that we don't find.
    /// @param scan_scopes indicares if we should continue the search for variabels in
    ///        successively outer (lower) scopes.
    /// @return An AST leaf node representing the variable found.
    [[nodiscard]] emp::Ptr<ASTNode_Var> ParseVar(ParseState & state,
                                                  bool create_ok=false,
                                                  bool scan_scopes=true);

    /// Load a value from the provided scope, which can come from a variable or a literal.
    [[nodiscard]] emp::Ptr<ASTNode> ParseValue(ParseState & state);

    /// Calculate the result of the provided operation on two computed entries.
    [[nodiscard]] emp::Ptr<ASTNode> ProcessOperation(const emp::Token & op_token,
                                                     emp::Ptr<ASTNode> value1,
                                                     emp::Ptr<ASTNode> value2);

    /// Calculate a full expression found in a token sequence, using the provided scope.
    /// @param state The current start of the parser and input stream
    /// @param decl_ok Can this expression begin with a declaration of a variable?
    /// @param prec_limit What is the highest precedence that expression should process?
    [[nodiscard]] emp::Ptr<ASTNode> ParseExpression(ParseState & state,
                                                    bool decl_ok=false,
                                                    size_t prec_limit=1000);

    /// Parse the declaration of a variable and return the newly created Symbol
    Var ParseDeclaration(ParseState & state);

    /// Parse an event description.
    emp::Ptr<ASTNode> ParseEvent(ParseState & state);

    /// Parse a specialty keyword statement (such as IF, WHILE, etc)
    emp::Ptr<ASTNode> ParseKeywordStatement(ParseState & state);

    /// Parse the next input in the specified Struct.  A statement can be a variable declaration,
    /// an expression, or an event.
    [[nodiscard]] emp::Ptr<ASTNode> ParseStatement(ParseState & state);

    /// Keep parsing statements until there aren't any more or we leave this scope. 
    [[nodiscard]] emp::Ptr<ASTNode_Block> ParseStatementList(ParseState & state) {
      Debug("Running ParseStatementList(", state.AsString(), ")");
      auto cur_block = emp::NewPtr<ASTNode_Block>(state.GetScope(), state.GetLine());
      cur_block->SetSymbolTable(state.GetSymbolTable());
      while (state.IsValid() && state.AsChar() != '}') {
        // Parse each statement in the file.
        emp::Ptr<ASTNode> statement_node = ParseStatement(state);

        // If the current statement is real, add it to the current block.
        if (!statement_node.IsNull()) cur_block->AddChild( statement_node );
      }
      return cur_block;
    }
  };

  // Load a variable name from the provided scope.
  emp::Ptr<ASTNode_Var> Parser::ParseVar(ParseState & state, bool create_ok, bool scan_scopes)
  {
    int start_line = state.GetLine();
    Debug("Running ParseVar(", state.AsString(), ",", create_ok, ",", scan_scopes,
          ") at line ", start_line);

    // TODO before the dots change we could use leading dots to search in a parent scope - is that something we need?

    // Next, we must have a variable name.
    // @CAO: Or a : ?  E.g., technically "..:size" could give you the parent scope size.
    state.RequireID("Must provide a variable identifier!");
    std::string var_name = state.UseLexeme();

    // Lookup this variable.
    Debug("...looking up symbol '", var_name,
          "' starting at scope '", state.GetScopeName(),
          "'; scanning=", scan_scopes);
    Var cur_symbol = state.LookupSymbol(var_name, scan_scopes);

    // If this variable just provided a scope, keep going.
    if (state.UseIfLexeme("::")) {
      state.PushScope(cur_symbol.GetValue()->AsScope());
      auto result = ParseVar(state, create_ok, false);
      state.PopScope();
      return result;
    }

    // Otherwise return the variable as a leaf!
    return emp::NewPtr<ASTNode_Var>(cur_symbol, var_name, start_line);
  }

  // Load a value from the provided scope, which can come from a variable or a literal.
  emp::Ptr<ASTNode> Parser::ParseValue(ParseState & state) {
    Debug("Running ParseValue(", state.AsString(), ")");

    // First check for a unary negation at the start of the value.
    if (state.UseIfChar('-')) {
      auto out_val = emp::NewPtr<ASTNode_Op1>("unary negation", state.GetLine());
      out_val->SetFun( [](double val){ return -val; } );
      out_val->AddChild(ParseValue(state));
      return out_val;
    }

    // Anything that begins with an identifier must represent a variable.  Refer!
    if (state.IsID()) return ParseVar(state, false, true);

    // A literal number should have a temporary created with its value.
    if (state.IsNumber()) {
      Debug("...value is a number: ", state.AsLexeme());
      double value = emp::from_string<double>(state.UseLexeme()); // Calculate value.
      return MakeTempLeaf(value);                                 // Return temporary Symbol.
    }

    // A literal string should be converted to a regular string and used.
    if (state.IsString()) {
      Debug("...value is a string: ", state.AsLexeme());
      std::string str = emp::from_literal_string(state.UseLexeme(), "\"'`"); // Convert literal string.
      return MakeTempLeaf(str);                                              // Return temporary Symbol.
    }

    // If we have an open parenthesis, process everything inside into a single value...
    if (state.UseIfChar('(')) {
      emp::Ptr<ASTNode> out_ast = ParseExpression(state);
      state.UseRequiredChar(')', "Expected a close parenthesis in expression.");
      return out_ast;
    }

    // Parse a list initializer [1, 2, 3]
    if (state.UseIfChar('[')) {
      auto list = emp::NewPtr<ASTNode_ListInit>(state.GetLine());
      while (!state.UseIfChar(']')) {
        list->AddChild(ParseExpression(state));
        if (!state.UseIfChar(',')) {
          state.UseRequiredChar(']', "Expected comma or ] after list element");
          break;
        }
      }
      return list;
    }

    state.Error("Expected a value, found: ", state.AsLexeme());

    return nullptr;
  }

  // Process a single provided operation on two Symbol objects.
  emp::Ptr<ASTNode> Parser::ProcessOperation(const emp::Token & op_token,
                                             emp::Ptr<ASTNode> in_node1,
                                             emp::Ptr<ASTNode> in_node2)
  {
    const std::string symbol = op_token.lexeme;
    emp_assert(!in_node1.IsNull());
    emp_assert(!in_node2.IsNull());

    // If this operation is assignment, do so!
    if (symbol == "=") return emp::NewPtr<ASTNode_Assign>(in_node1, in_node2, op_token.line_id);
    
    // Determine the output value and put it in a temporary node.
    emp::Ptr<ASTNode_Op2> out_val = emp::NewPtr<ASTNode_Op2>(symbol, op_token.line_id);

    if (symbol == "+") out_val->SetFun( [](emp::Datum v1, emp::Datum v2){ return v1 + v2; } );
    else if (symbol == "-") out_val->SetFun( [](emp::Datum v1, emp::Datum v2){ return v1 - v2; } );
    else if (symbol == "**") out_val->SetFun( [](emp::Datum v1, emp::Datum v2){ return emp::Pow(v1.AsDouble(), v2.AsDouble()); } );
    else if (symbol == "*") out_val->SetFun( [](emp::Datum v1, emp::Datum v2){ return v1 * v2; } );
    else if (symbol == "/") out_val->SetFun( [](emp::Datum v1, emp::Datum v2){ return v1 / v2; } );
    else if (symbol == "%") out_val->SetFun( [](emp::Datum v1, emp::Datum v2){ return v1 % v2; } );
    else if (symbol == "==") out_val->SetFun( [](emp::Datum v1, emp::Datum v2){ return v1 == v2; } );
    else if (symbol == "!=") out_val->SetFun( [](emp::Datum v1, emp::Datum v2){ return v1 != v2; } );
    else if (symbol == "<")  out_val->SetFun( [](emp::Datum v1, emp::Datum v2){ return v1 < v2; } );
    else if (symbol == "<=") out_val->SetFun( [](emp::Datum v1, emp::Datum v2){ return v1 <= v2; } );
    else if (symbol == ">")  out_val->SetFun( [](emp::Datum v1, emp::Datum v2){ return v1 > v2; } );
    else if (symbol == ">=") out_val->SetFun( [](emp::Datum v1, emp::Datum v2){ return v1 >= v2; } );

    // @CAO: Need to still handle these last two differently for short-circuiting.
    else if (symbol == "&&") out_val->SetFun( [](emp::Datum v1, emp::Datum v2){ return v1 && v2; } );
    else if (symbol == "||") out_val->SetFun( [](emp::Datum v1, emp::Datum v2){ return v1 || v2; } );

    out_val->AddChild(in_node1);
    out_val->AddChild(in_node2);

    return out_val;
  }
                                      

  /// Calculate a full expression found in a token sequence, using the provided scope.
  /// @param state The current start of the parser and input stream
  /// @param decl_ok Can this expression begin with a declaration of a variable?
  /// @param prec_limit What is the highest precedence that expression should process?

  emp::Ptr<ASTNode> Parser::ParseExpression(ParseState & state, bool decl_ok, size_t prec_limit) {
    Debug("Running ParseExpression(", state.AsString(), ", decl_ok=", decl_ok, ", limit=", prec_limit, ")");

    // Allow this statement to be a declaration if it begins with a type.
    if (decl_ok && state.IsType()) {
      Var new_symbol = ParseDeclaration(state);

      // Functions can't be used while they're being defined, so return them now
      if (new_symbol.GetValue()->IsFunction()) {
        return emp::NewPtr<ASTNode_Var>(new_symbol, new_symbol.GetValue()->GetName());
      }
 
      // If this symbol is a new scope, it can be populated now either directly (in braces)
      // or indirectly (with an assignment)
      if (new_symbol.GetValue()->IsScope()) {
        if (state.UseIfChar('{')) {
          // Make a new scope with an inaccessible name to use for initialization
          emp::Ptr<Symbol_Scope> scope = state
            .AddScope(new_symbol.GetValue()->GetName() + "'", "Local struct initialization scope")
            .GetValue()->AsScopePtr();
          state.PushScope(*scope);
          emp::Ptr<ASTNode_Block> out_node = ParseStatementList(state);
          state.PopScope();

          // At the end of the block, assign the variable to a copy of the created scope
          // This way, running the same struct declaration twice results in separate struct instances
          auto clone_node = emp::NewPtr<ASTNode_Clone>();
          clone_node->AddChild(emp::NewPtr<ASTNode_Leaf>(scope));
          auto var_node = emp::NewPtr<ASTNode_Var>(new_symbol, scope->GetName());
          out_node->AddChild(emp::NewPtr<ASTNode_Assign>(var_node, clone_node));

          state.UseRequiredChar('}', "Expected scope '", scope->GetName(), "' to end with a '}'.");
          return out_node;
        }
      }      

      // Otherwise rewind so that the new variable can be used to start an expression.
      --state;
    }

    /// Process a value (and possibly more!)
    emp::Ptr<ASTNode> cur_node = ParseValue(state);

    while (state.UseIfChar('.')) {
      state.RequireID("Expected member name after '.'");
      std::string name = state.UseLexeme();
      auto node = emp::NewPtr<ASTNode_Member>(name);
      node->AddChild(cur_node);
      cur_node = node;
    }

    emp::Token op_token = state.AsToken();
    std::string op = state.AsLexeme();

    Debug("...back in ParseExpression; op=`", op, "`; state=", state.AsString());

    while ( emp::Has(precedence_map, op) && precedence_map[op] < prec_limit ) {
      ++state; // Move past the current operator
      // Do we have a function call?
      if (op == "(") {
        // Collect arguments.
        emp::vector< emp::Ptr<ASTNode> > args;
        while (state.AsChar() != ')') {
          emp::Ptr<ASTNode> next_arg = ParseExpression(state);
          args.push_back(next_arg);          // Save this argument.
          if (state.AsChar() != ',') break;  // If we don't have a comma, no more args!
          ++state;                           // Move on to the next argument.
        }
        state.UseRequiredChar(')', "Expected a ')' to end function call.");

        // cur_node should have evaluated itself to a function; a Call node will link that
        // function with its arguments, run it, and return the result.
        cur_node = emp::NewPtr<ASTNode_Call>(cur_node, args, op_token.line_id);
      }
      // Do we have an array subscript?
      else if (op == "[") {
        auto idx_node = ParseExpression(state);
        state.UseRequiredChar(']', "Expected a ']' to end subscript index.");
        auto node = emp::NewPtr<ASTNode_Subscript>();
        node->AddChild(cur_node);
        node->AddChild(idx_node);
        cur_node = node;
      }
      // Otherwise we must have a binary math operation.
      else {
        emp::Ptr<ASTNode> node2 = ParseExpression(state, false, precedence_map[op]);
        cur_node = ProcessOperation(op_token, cur_node, node2);
      }

      // Move the current value over to cur_node and check if we have a new operator...
      op = state.AsLexeme();
      op_token = state.AsToken();
    }

    emp_assert(!cur_node.IsNull());
    return cur_node;
  }

  // Parse an the declaration of a variable.
  Var Parser::ParseDeclaration(ParseState & state) {
    std::string type_name = state.UseLexeme();
    state.RequireID("Type name '", type_name, "' must be followed by variable to declare.");
    std::string var_name = state.UseLexeme();

    if (state.AsLexeme() == "(") {
      // This is a function
      state.UseLexeme();
      // Add a space to the scope name so it can't be used from the script
      emp::Ptr<Symbol_Scope> scope = emp::NewPtr<Symbol_Scope>(var_name, "Function scope", &state.GetScope());
      state.PushScope(*scope);

      emp::vector<Var> params;
      while (state.AsChar() != ')') {
        state.RequireID("Function parameters must be names (without types).");
        std::string param_name = state.UseLexeme();
        Var param = state.AddLocalVar(param_name, "Function parameter.");
        params.push_back(param);
        state.UseIfChar(',');                     // Skip comma if next (does allow trailing comma)
      }
      state.UseRequiredChar(')', "Function args must end in a ')'");

      state.UseRequiredChar('{', "Function body must start with '{'");
      auto body_block = ParseStatementList(state);
      state.UseRequiredChar('}', "Function body must end with '{'");
      state.PopScope();

      return state.AddFunction(var_name, "Local function.", type_name, params, body_block, scope);
    }

    if (type_name == "Var") return state.AddLocalVar(var_name, "Local variable.");
    else if (type_name == "Struct") return state.AddScope(var_name, "Local struct");

    // Otherwise we have an object of a custom type to add.
    Debug("Building object '", var_name, "' of type '", type_name, "'");
    return state.AddObject(type_name, var_name);
  }

  // Parse an event description.
  emp::Ptr<ASTNode> Parser::ParseEvent(ParseState & state) {
    emp::Token start_token = state.AsToken();
    state.UseRequiredChar('@', "All event declarations must being with an '@'.");
    state.RequireID("Events must start by specifying signal name.");
    const std::string & trigger_name = state.UseLexeme();
    state.UseRequiredChar('(', "Expected parentheses after '", trigger_name, "' for args.");

    emp::vector<emp::Ptr<ASTNode>> args;
    while (state.AsChar() != ')') {
      args.push_back( ParseExpression(state, true) );
      state.UseIfChar(',');                     // Skip comma if next (does allow trailing comma)
    }
    state.UseRequiredChar(')', "Event args must end in a ')'");

    auto action_block = emp::NewPtr<ASTNode_Block>(state.GetScope(), state.GetLine());
    action_block->SetSymbolTable(state.GetSymbolTable());
    emp::Ptr<ASTNode> action_node = ParseStatement(state);

    // If the action statement is real, add it to the action block.
    if (!action_node.IsNull()) action_block->AddChild( action_node );

    Debug("Building event '", trigger_name, "' with args ", args);

    state.AddAction(trigger_name, args, action_block, start_token.line_id);

    return nullptr;
  }

  /// Parse a specialty keyword statement (such as IF, WHILE, etc)
  emp::Ptr<ASTNode> Parser::ParseKeywordStatement(ParseState & state) {
    size_t keyword_line = state.GetLine();

    if (state.UseIfLexeme("IF")) {
      state.UseRequiredChar('(', "Expected '(' to begin IF test condition.");
      emp::Ptr<ASTNode> test_node = ParseExpression(state);
      state.UseRequiredChar(')', "Expected ')' to end IF test condition.");
      emp::Ptr<ASTNode> true_node = ParseStatement(state);
      emp::Ptr<ASTNode> else_node = nullptr;
      if (state.UseIfLexeme("ELSE")) else_node = ParseStatement(state);
      return emp::NewPtr<ASTNode_If>(test_node, true_node, else_node, keyword_line);
    }

    else if (state.UseIfLexeme("WHILE")) {
      state.UseRequiredChar('(', "Expected '(' to begin IF test condition.");
      emp::Ptr<ASTNode> test_node = ParseExpression(state);
      state.UseRequiredChar(')', "Expected ')' to end IF test condition.");
      emp::Ptr<ASTNode> body_node = ParseStatement(state);
      return emp::NewPtr<ASTNode_While>(test_node, body_node, keyword_line);
    }

    else if (state.UseIfLexeme("BREAK")) { return MakeBreakLeaf(keyword_line); }

    else if (state.UseIfLexeme("CONTINUE")) { return MakeContinueLeaf(keyword_line); }

    else if (state.UseIfLexeme("RETURN")) {
      auto node = emp::NewPtr<ASTNode_Return>(keyword_line);
      // Allow return without a value if it's followed by a semicolon
      if (state.AsChar() != ';')
        node->AddChild(ParseExpression(state));
      return node;
    }

    // If we made it this far, we have an error.  Identify and deal with it!

    if (state.UseIfLexeme("ELSE")) state.Error("'ELSE' must be preceded by an 'IF' statement.");
    else state.Error("Keyword '", state.AsLexeme(), "' not yet implemented.");

    return nullptr;
  }

  // Process the next input in the specified Struct.
  emp::Ptr<ASTNode> Parser::ParseStatement(ParseState & state) {
    Debug("Running ParseStatement(", state.AsString(), ")");

    // Allow a statement with an empty line.
    if (state.UseIfChar(';')) { return nullptr; }

    // Allow a statement to be a new scope.
    if (state.UseIfChar('{')) {
      // @CAO Need to add an anonymous scope (that gets written properly)
      emp::Ptr<ASTNode> out_node = ParseStatementList(state);
      state.UseRequiredChar('}', "Expected '}' to close scope.");
      return out_node;
    }

    // Allow event definitions if a statement begins with an '@'
    if (state.AsChar() == '@') return ParseEvent(state);

    // Allow select commands that are only possible at the full statement level (not expressions)
    if (state.IsKeyword()) return ParseKeywordStatement(state);

    // If we made it here, remainder should be an expression; it may begin with a declaration.
    emp::Ptr<ASTNode> out_node = ParseExpression(state, true);

    // Expressions must end in a semi-colon.
    state.UseRequiredChar(';', "Expected ';' at the end of a statement; found: ", state.AsLexeme());

    return out_node;
  }  
}
#endif
