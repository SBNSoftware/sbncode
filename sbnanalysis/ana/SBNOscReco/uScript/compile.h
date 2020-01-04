#ifndef uscript_compile_h
#define uscript_compile_h

#include <unordered_set>

#include "vm.h"
#include "chunk.h"
#include "scanner.h"
#include "tclass.h"

namespace uscript {

class Compiler {
  struct Parser {
    Token current;
    Token previous;
    bool hadError;
    bool panicMode;
  };

  enum Precedence {
    PREC_NONE,
    PREC_ASSIGNMENT,  // =
    PREC_OR,          // or
    PREC_AND,         // and
    PREC_EQUALITY,    // == !=
    PREC_COMPARISON,  // < > <= >=
    PREC_TERM,        // + -
    PREC_FACTOR,      // * /
    PREC_UNARY,       // ! -
    PREC_CALL,        // . ()
    PREC_PRIMARY
  };

  using ParseFn = void (Compiler::*)(bool);
  // typedef void (Compiler::*ParseFn)();

  struct ParseRule {
    ParseFn prefix;
    ParseFn infix;
    Precedence precedence;
  };

  struct Local {
    Token name;
    int depth;
  };

  std::unordered_set<std::string> strings;
  Parser parser;
  Scanner scanner;
  Chunk *current;
  std::vector<Local> locals;
  int scopeDepth;

  TClassList tclasslist;

  void Declaration();
  void Statement();
  void Synchronize();
  void Block();

  void BeginScope();
  void EndScope();

  void VarDeclaration();
  uint8_t IdentifierConstant(Token *token);
  uint8_t ParseVariable(const char *errorMessage);
  void DefineVariable(uint8_t global);
  void DeclareVariable();
  void AddLocal(Token name);
  int ResolveLocal(Token *name);
  void MarkIntialized();
  void PrintStatement();
  void ExpressionStatement();
  void IfStatement();
  void WhileStatement();
  void ForStatement();

  bool Match(TokenType type);
  bool Check(TokenType type);

  void ErrorAt(Token &token, const char *message);
  void PatchJump(int offset);
  int EmitJump(uint8_t instruction);
  void EmitLoop(int loopStart);
  void EmitByte(uint8_t byte);
  void EmitBytes(uint8_t bytea, uint8_t byteb);
  void EmitReturn();
  void EmitConstant(Value value);

  uint8_t MakeConstant(Value value);

  void Number(bool canAssign);
  void Grouping(bool canAssign);
  void Unary(bool canAssign);
  void Binary(bool canAssign);
  void Literal(bool canAssign);
  void String(bool canAssign);
  void Variable(bool canAssign);
  void And(bool canAssign);
  void Or(bool canAssign);
  void Call(bool canAssign);
  void Dot(bool canAssign);

  void NamedVariable(Token token, bool canAssign);

  void ParsePrecedence(Precedence precedence);

  ParseRule *GetRule(TokenType type) const;

  uint8_t ArgumentList();

  void Advance();
  void Expression();
  void Consume(TokenType type, const char *message);
  void SetChunk(Chunk *chunk) { current = chunk; }
  void Finish();
  Chunk *CurrentChunk() { return current; }

  Compiler();

  bool DoCompile(const char *source, Chunk *chunk);
  const char *DoIntern(const std::string &str);
  void DoRegister(const char *classname);

  static Compiler &Instance() {
    static Compiler c;
    return c;
  }

public:
  static const char *Intern(const std::string &str) { return Instance().DoIntern(str); }
  static bool HadError() { return Instance().parser.hadError; }
  static bool Compile(const char *source, Chunk *chunk) { Instance().DoCompile(source, chunk); }
  
  static TClassInfo *GetClassInfo(const char *classname) {
    return Instance().tclasslist.Add(classname);
  }

  template <typename TObj>
  static void inline Register() {
    Instance().DoRegister(std::string(type_name<TObj>()).c_str());
  }

};


} // end namespace
#endif
