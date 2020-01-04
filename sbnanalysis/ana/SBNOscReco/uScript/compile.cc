#include <iostream>
#include <string.h>
#include <cassert>

#include "compile.h"
#include "scanner.h"

bool uscript::Compiler::DoCompile(const char *source, Chunk *chunk) {
  scanner.SetSource(source);

  SetChunk(chunk);

  Advance();

  while (!Match(uscript::TOKEN_EOF)) {
    Declaration();
  }

  Consume(uscript::TOKEN_EOF, "Expect end of expression.");

  Finish();

  return !parser.hadError;
}

uscript::Compiler::Compiler() {
  parser.hadError = false;
  parser.panicMode = false;
  scopeDepth = 0;
}

void uscript::Compiler::Declaration() {
  if (Match(uscript::TOKEN_VAR)) {
    VarDeclaration();
  }
  else {
    Statement();
  }

  if (parser.panicMode) Synchronize();
}

void uscript::Compiler::VarDeclaration() {
  uint8_t global = ParseVariable("Expect variable name.");

  if (Match(uscript::TOKEN_EQUAL)) {
    Expression();
  }
  else {
    EmitByte(uscript::OP_NIL);
  }

  Consume(uscript::TOKEN_SEMICOLON, "Expect ';' after declaration");

  DefineVariable(global);
}

uint8_t uscript::Compiler::IdentifierConstant(uscript::Token *token) {
  std::string str(token->start, token->start + token->length);
  const char *interned = DoIntern(str);
  return MakeConstant(STRING_VAL(interned));
}

uint8_t uscript::Compiler::ParseVariable(const char *errorMessage) {
  Consume(uscript::TOKEN_IDENTIFIER, errorMessage);
  DeclareVariable();
  if (scopeDepth > 0) return 0; // local variable 
  return IdentifierConstant(&parser.previous); // global variable
}

static bool identifiersEqual(uscript::Token* a, uscript::Token* b) {
  if (a->length != b->length) return false;
  return memcmp(a->start, b->start, a->length) == 0;
}

int uscript::Compiler::ResolveLocal(Token *name) {
  for (int i = locals.size() -1; i >= 0; i--) {
    if (identifiersEqual(name, &locals[i].name)) {
      if (locals[i].depth == -1) {
        ErrorAt(parser.previous, "Cannot read from local variable in its own initializer.");
      }
      return i;
    }
  } 
  return -1;
}

void uscript::Compiler::DeclareVariable() {
  if (scopeDepth == 0) return; // global variables are implicitly declared

  Token name = parser.previous;
  // make sure the local isn't already there
  for (int i = locals.size()-1; i >= 0; i--) {
    if (locals[i].depth != -1 && locals[i].depth < scopeDepth) {
      break;
    }
    if (identifiersEqual(&locals[i].name, &name)) {
      ErrorAt(parser.previous, "Variable with this name already declared in this scope.");
    }
  }

  AddLocal(name);
}

void uscript::Compiler::AddLocal(Token name) {
  if (locals.size() > UINT8_MAX) {
    ErrorAt(parser.previous, "Too many local variables in one chunk.");
  }
  Local local;
  local.name = name;
  local.depth = -1;
  locals.push_back(local);
}

void uscript::Compiler::MarkIntialized() {
  locals[locals.size() -1].depth = scopeDepth;
}

void uscript::Compiler::DefineVariable(uint8_t global) {
  if (scopeDepth > 0) {
    MarkIntialized();
    return;
  }
  EmitBytes(uscript::OP_DEFINE_GLOBAL, global);
}


void uscript::Compiler::Synchronize() {
  parser.panicMode = false;

  while (parser.current.type != TOKEN_EOF) {
    if (parser.previous.type == TOKEN_SEMICOLON) return;

    switch (parser.current.type) {
      case TOKEN_CLASS:
      case TOKEN_FUN:
      case TOKEN_VAR:
      case TOKEN_FOR:
      case TOKEN_IF:
      case TOKEN_WHILE:
      case TOKEN_PRINT:
      case TOKEN_RETURN:
        return;
      default:
        break;
    }
    Advance();
  }
}

void uscript::Compiler::Statement() {
  if (Match(uscript::TOKEN_PRINT)) {
    PrintStatement();
  }
  else if (Match(uscript::TOKEN_IF)) {
    IfStatement();
  }
  else if (Match(uscript::TOKEN_WHILE)) {
    WhileStatement();
  }
  else if (Match(uscript::TOKEN_FOR)) {
    ForStatement();
  }
  else if (Match(uscript::TOKEN_LEFT_BRACE)) {
    BeginScope();
    Block();
    EndScope();
  }
  else {
    ExpressionStatement();
  }
}

void uscript::Compiler::WhileStatement() {
  int loopStart = CurrentChunk()->code.size();

  Consume(uscript::TOKEN_LEFT_PAREN, "Expect '(' after while.");
  Expression();
  Consume(uscript::TOKEN_RIGHT_PAREN, "Expect ')' after condition.");

  int exitJump = EmitJump(uscript::OP_JUMP_IF_FALSE);

  EmitByte(uscript::OP_POP);
  Statement();

  PatchJump(exitJump);
  EmitByte(uscript::OP_POP);

  EmitLoop(loopStart);
}

void uscript::Compiler::ForStatement() {
  BeginScope();

  Consume(uscript::TOKEN_LEFT_PAREN, "Expect '(' after 'for'.");
  if (Match(uscript::TOKEN_SEMICOLON)) {
    // no initializer
  }
  else if (Match(uscript::TOKEN_VAR)) {
    VarDeclaration();
  }
  else {
    ExpressionStatement();
  }

  int loopStart = CurrentChunk()->code.size();

  int exitJump = -1;
  if (!Match(uscript::TOKEN_SEMICOLON)) {
    Expression();
    Consume(uscript::TOKEN_SEMICOLON, "Expect ';' after loop condition.");
    exitJump = EmitJump(uscript::OP_JUMP_IF_FALSE);
    EmitByte(uscript::OP_POP);
  }

  if (!Match(uscript::TOKEN_RIGHT_PAREN)) {
    int bodyJump = EmitJump(uscript::OP_JUMP);
    int incrementStart = CurrentChunk()->code.size();

    Expression();
    EmitByte(uscript::OP_POP);

    Consume(uscript::TOKEN_RIGHT_PAREN, "Expect ')' after for clause.");

    EmitLoop(loopStart);
    loopStart = incrementStart;
    PatchJump(bodyJump);
    
  }

  // body of for loop
  Statement();

  EmitLoop(loopStart);

  if (exitJump != -1) {
    PatchJump(exitJump);
    EmitByte(uscript::OP_POP);
  }
  
  EndScope();
}

void uscript::Compiler::IfStatement() {
  Consume(uscript::TOKEN_LEFT_PAREN, "Expect '(' after if.");
  Expression();
  Consume(uscript::TOKEN_RIGHT_PAREN, "Expect ')' after condition.");

  int thenJump = EmitJump(uscript::OP_JUMP_IF_FALSE);
  EmitByte(uscript::OP_POP);
  Statement();
  PatchJump(thenJump);

  int elseJump = EmitJump(uscript::OP_JUMP);
  if (Match(uscript::TOKEN_ELSE)) Statement();
  PatchJump(elseJump);
}

int uscript::Compiler::EmitJump(uint8_t instruction) {
  EmitByte(instruction);
  EmitByte(0xff);
  EmitByte(0xff);
  return CurrentChunk()->code.size() - 2;
}

void uscript::Compiler::PatchJump(int offset) {
  int jump = CurrentChunk()->code.size() - offset - 2;

  if (jump > UINT16_MAX) {
    ErrorAt(parser.previous, "Too much code to jump over.");
  }

  CurrentChunk()->code[offset] = (jump >> 8) & 0xff; 
  CurrentChunk()->code[offset + 1] = jump & 0xff;
}

void uscript::Compiler::Block() {
  while (!Check(uscript::TOKEN_RIGHT_BRACE) && !Check(uscript::TOKEN_EOF)) {
    Declaration();
  }
  Consume(uscript::TOKEN_RIGHT_BRACE, "Expect '}' after block.");
}

void uscript::Compiler::BeginScope() {
  scopeDepth ++;
}

void uscript::Compiler::EndScope() {
  scopeDepth --;

  while (locals.size() > 0 && locals[locals.size()-1].depth > scopeDepth) {
    EmitByte(uscript::OP_POP);
    locals.pop_back();
  } 

}

void uscript::Compiler::PrintStatement() {
  Expression();
  Consume(uscript::TOKEN_SEMICOLON, "Expect ';' after value.");
  EmitByte(uscript::OP_PRINT);
}

void uscript::Compiler::ExpressionStatement() {
  Expression();
  Consume(uscript::TOKEN_SEMICOLON, "Expect ';' after value.");
  EmitByte(uscript::OP_POP);
}

bool uscript::Compiler::Match(uscript::TokenType type) {
  if (!Check(type)) return false;
  Advance();
  return true;
}

bool uscript::Compiler::Check(uscript::TokenType type) {
  return parser.current.type == type;
}

void uscript::Compiler::EmitLoop(int loopStart) {
  EmitByte(uscript::OP_LOOP);

  int offset = CurrentChunk()->code.size() - loopStart + 2;

  if (offset > UINT16_MAX) ErrorAt(parser.previous, "Loop body too large.");

  EmitByte((offset >> 8) & 0xff);
  EmitByte(offset & 0xff);
}

void uscript::Compiler::EmitByte(uint8_t byte) {
  CurrentChunk()->Write(byte);
}

void uscript::Compiler::EmitBytes(uint8_t bytea, uint8_t byteb) {
  CurrentChunk()->Write(bytea);
  CurrentChunk()->Write(byteb);
}

static void PrintErrorAt(uscript::Token &token, const char *message) {
  std::cerr << "Error";
  if (token.type == uscript::TOKEN_EOF) {
    std::cerr << " at end";
  }
  else if (token.type == uscript::TOKEN_ERROR) {
    // nothing 
  }
  else {
    std::string error(token.start, token.start + token.length);
    std::cerr << " at " << error ;
  }

  std::cerr << ": " << message << std::endl;
  return;
}

void uscript::Compiler::ErrorAt(uscript::Token &token, const char *message) {
  if (parser.panicMode) return;

  parser.hadError = true;
  parser.panicMode = true;
  PrintErrorAt(token, message);
}

void uscript::Compiler::Expression() {
  ParsePrecedence(uscript::Compiler::PREC_ASSIGNMENT);
}

void uscript::Compiler::Advance() {
  parser.previous = parser.current;

  while (1) {
    parser.current = scanner.ScanToken();
    if (parser.current.type != uscript::TOKEN_ERROR) break;

    ErrorAt(parser.current, parser.current.start);
  }
}

void uscript::Compiler::Consume(uscript::TokenType type, const char *message) {
  if (parser.current.type == type) {
    Advance();
    return;
  }

  ErrorAt(parser.current, message);
}

void uscript::Compiler::Finish() {
  EmitReturn();
#ifdef DEBUG_PRINT_CODE
  if (!parser.hadError) {
    CurrentChunk()->Disassemble("code");
  }
#endif
}

void uscript::Compiler::EmitReturn() {
  return EmitByte(uscript::OP_RETURN);
}

void uscript::Compiler::EmitConstant(Value value) {
  return EmitBytes(uscript::OP_CONSTANT, MakeConstant(value));
}

uint8_t uscript::Compiler::MakeConstant(Value value) {
  int constant = CurrentChunk()->AddConstant(value);
  if (constant > UINT8_MAX) {
    ErrorAt(parser.previous, "Too many constants in one chunk.");
    return 0;
  }
  return constant;
}

void uscript::Compiler::Number(bool canAssign) {
  // try to see if we can convert to integer
  char *end;
  int v_int = strtol(parser.previous.start, &end, 10 /* base 10*/);
  if (end == parser.previous.start + parser.previous.length) {
    return EmitConstant(INTEGER_VAL(v_int));
  }
  else {
    double value = strtod(parser.previous.start, NULL);
    EmitConstant(NUMBER_VAL(value));
  }
}

void uscript::Compiler::Grouping(bool canAssign) {
  Expression();
  Consume(uscript::TOKEN_RIGHT_PAREN, "Expect ')' after expression.");
}

void uscript::Compiler::And(bool canAssign) {
  int endJump = EmitJump(uscript::OP_JUMP_IF_FALSE);

  EmitByte(uscript::OP_POP);

  ParsePrecedence(PREC_AND);

  PatchJump(endJump);
}

uint8_t uscript::Compiler::ArgumentList() {
  uint8_t argCount = 0;
  if (!Check(uscript::TOKEN_RIGHT_PAREN)) {
    do {
      Expression();
      if (argCount == 255) {
        ErrorAt(parser.previous, "Cannot have more than 255 arguments.");
      }
      argCount++;
    } while (Match(uscript::TOKEN_COMMA));
  }

  Consume(uscript::TOKEN_RIGHT_PAREN, "Expect ')' after function arguments.");
  return argCount;
}

void uscript::Compiler::Dot(bool canAssign) {
  Consume(uscript::TOKEN_IDENTIFIER, "Expect property name after '.'.");
  uint8_t name = IdentifierConstant(&parser.previous); 

  EmitBytes(uscript::OP_GET_PROPERTY, name);
}

void uscript::Compiler::Call(bool canAssign) {
  uint8_t argCount = ArgumentList();
  EmitBytes(uscript::OP_CALL, argCount);
}

void uscript::Compiler::Or(bool canAssign) {
  int elseJump = EmitJump(uscript::OP_JUMP_IF_FALSE);
  int endJump =  EmitJump(uscript::OP_JUMP);

  PatchJump(elseJump);
  EmitByte(uscript::OP_POP);

  ParsePrecedence(PREC_OR);
  PatchJump(endJump);
}

void uscript::Compiler::Unary(bool canAssign) {
  uscript::TokenType operatorType = parser.previous.type;

  // compile the operand
  ParsePrecedence(uscript::Compiler::PREC_UNARY);

  switch (operatorType) {
    case uscript::TOKEN_BANG: EmitByte(uscript::OP_NOT); break;
    case uscript::TOKEN_MINUS: EmitByte(uscript::OP_NEGATE); break;
    default: return; // unreachable
  }
}

void uscript::Compiler::Literal(bool canAssign) {
  switch (parser.previous.type) {
    case uscript::TOKEN_FALSE: EmitByte(uscript::OP_FALSE); break;
    case uscript::TOKEN_TRUE: EmitByte(uscript::OP_TRUE); break;
    case uscript::TOKEN_NIL: EmitByte(uscript::OP_NIL); break;
    default:
      return; //unreachable
  }
}

void uscript::Compiler::String(bool canAssign) {
  std::string str(parser.previous.start+1, parser.previous.start + parser.previous.length - 1);
  EmitConstant(STRING_VAL(DoIntern(str)));
}

void uscript::Compiler::Variable(bool canAssign) {
  NamedVariable(parser.previous, canAssign);
}

void uscript::Compiler::NamedVariable(uscript::Token name, bool canAssign) {
  uint8_t getOp, setOp;
  int arg = ResolveLocal(&name);
  if (arg != -1) {
    getOp = OP_GET_LOCAL;
    setOp = OP_SET_LOCAL;
  } 
  else {
    arg = IdentifierConstant(&name);
    getOp = OP_GET_GLOBAL;
    setOp = OP_SET_GLOBAL;
  }

  if (canAssign && Match(uscript::TOKEN_EQUAL)) {
    Expression();
    EmitBytes(setOp, (uint8_t)arg);
  }
  else {
    EmitBytes(getOp, (uint8_t)arg);
  }
}

void uscript::Compiler::Binary(bool canAssign) {
  uscript::TokenType operatorType = parser.previous.type;

  uscript::Compiler::ParseRule *rule = GetRule(operatorType);

  ParsePrecedence((uscript::Compiler::Precedence) (rule->precedence + 1));

  switch (operatorType) {
    case uscript::TOKEN_PLUS:    EmitByte(uscript::OP_ADD); break;
    case uscript::TOKEN_MINUS:   EmitByte(uscript::OP_SUBTRACT); break;
    case uscript::TOKEN_STAR:    EmitByte(uscript::OP_MULTIPLY); break;
    case uscript::TOKEN_SLASH:   EmitByte(uscript::OP_DIVIDE); break;
    case uscript::TOKEN_BANG_EQUAL:    EmitBytes(uscript::OP_EQUAL, uscript::OP_NOT); break;
    case uscript::TOKEN_EQUAL_EQUAL:   EmitByte(uscript::OP_EQUAL); break;
    case uscript::TOKEN_GREATER:       EmitByte(uscript::OP_GREATER); break;
    case uscript::TOKEN_GREATER_EQUAL: EmitBytes(uscript::OP_LESS, uscript::OP_NOT); break;
    case uscript::TOKEN_LESS:          EmitByte(uscript::OP_LESS); break;
    case uscript::TOKEN_LESS_EQUAL:    EmitBytes(uscript::OP_GREATER, uscript::OP_NOT); break;
    default:
      return; // unreachable
  }
}

void uscript::Compiler::ParsePrecedence(uscript::Compiler::Precedence precedence) {
  Advance();

  ParseFn prefixRule = GetRule(parser.previous.type)->prefix;
  if (prefixRule == NULL) {
    ErrorAt(parser.previous, "Expect expression.");
    return;
  }
  
  bool canAssign = precedence <= uscript::Compiler::PREC_ASSIGNMENT;
  (*this.*prefixRule)(canAssign);

  while (precedence <= GetRule(parser.current.type)->precedence) {
    Advance();
    ParseFn infixRule = GetRule(parser.previous.type)->infix;
    (*this.*infixRule)(canAssign);
  }

  if (canAssign && Match(uscript::TOKEN_EQUAL)) {
    ErrorAt(parser.previous, "Invalid assignmnet target.");
  }
}

uscript::Compiler::ParseRule *uscript::Compiler::GetRule(uscript::TokenType type) const {
  static ParseRule rules[] = {
    { &uscript::Compiler::Grouping, &uscript::Compiler::Call,    uscript::Compiler::PREC_NONE },       // TOKEN_LEFT_PAREN
    { NULL,     NULL,    uscript::Compiler::PREC_NONE },       // TOKEN_RIGHT_PAREN
    { NULL,     NULL,    uscript::Compiler::PREC_NONE },       // TOKEN_LEFT_BRACE
    { NULL,     NULL,    uscript::Compiler::PREC_NONE },       // TOKEN_RIGHT_BRACE
    { NULL,     NULL,    uscript::Compiler::PREC_NONE },       // TOKEN_COMMA
    { NULL,     &uscript::Compiler::Dot,    uscript::Compiler::PREC_CALL },       // TOKEN_DOT
    { &uscript::Compiler::Unary,    &uscript::Compiler::Binary,  uscript::Compiler::PREC_TERM },       // TOKEN_MINUS
    { NULL,     &uscript::Compiler::Binary,  uscript::Compiler::PREC_TERM },       // TOKEN_PLUS
    { NULL,     NULL,    uscript::Compiler::PREC_NONE },       // TOKEN_SEMICOLON
    { NULL,     &uscript::Compiler::Binary,  uscript::Compiler::PREC_FACTOR },     // TOKEN_SLASH
    { NULL,     &uscript::Compiler::Binary,  uscript::Compiler::PREC_FACTOR },     // TOKEN_STAR
    { &uscript::Compiler::Unary,     NULL,    uscript::Compiler::PREC_NONE },       // TOKEN_BANG
    { NULL,     &uscript::Compiler::Binary,    uscript::Compiler::PREC_EQUALITY },       // TOKEN_BANG_EQUAL
    { NULL,     NULL,    uscript::Compiler::PREC_NONE },       // TOKEN_EQUAL
    { NULL,     &uscript::Compiler::Binary,    uscript::Compiler::PREC_EQUALITY },       // TOKEN_EQUAL_EQUAL
    { NULL,     &uscript::Compiler::Binary,    uscript::Compiler::PREC_COMPARISON },       // TOKEN_GREATER
    { NULL,     &uscript::Compiler::Binary,    uscript::Compiler::PREC_COMPARISON },       // TOKEN_GREATER_EQUAL
    { NULL,     &uscript::Compiler::Binary,    uscript::Compiler::PREC_COMPARISON },       // TOKEN_LESS
    { NULL,     &uscript::Compiler::Binary,    uscript::Compiler::PREC_COMPARISON },       // TOKEN_LESS_EQUAL
    { &uscript::Compiler::Variable,     NULL,    uscript::Compiler::PREC_NONE },       // TOKEN_IDENTIFIER
    { &uscript::Compiler::String,     NULL,    uscript::Compiler::PREC_NONE },       // TOKEN_STRING
    { &uscript::Compiler::Number,   NULL,    uscript::Compiler::PREC_NONE },       // TOKEN_NUMBER
    { NULL,     &uscript::Compiler::And,    uscript::Compiler::PREC_AND },       // TOKEN_AND
    { NULL,     NULL,    uscript::Compiler::PREC_NONE },       // TOKEN_CLASS
    { NULL,     NULL,    uscript::Compiler::PREC_NONE },       // TOKEN_ELSE
    { &uscript::Compiler::Literal,     NULL,    uscript::Compiler::PREC_NONE },       // TOKEN_FALSE
    { NULL,     NULL,    uscript::Compiler::PREC_NONE },       // TOKEN_FOR
    { NULL,     NULL,    uscript::Compiler::PREC_NONE },       // TOKEN_FUN
    { NULL,     NULL,    uscript::Compiler::PREC_NONE },       // TOKEN_IF
    { &uscript::Compiler::Literal,     NULL,    uscript::Compiler::PREC_NONE },       // TOKEN_NIL
    { NULL,     &uscript::Compiler::Or,    uscript::Compiler::PREC_OR },       // TOKEN_OR
    { NULL,     NULL,    uscript::Compiler::PREC_NONE },       // TOKEN_PRINT
    { NULL,     NULL,    uscript::Compiler::PREC_NONE },       // TOKEN_RETURN
    { NULL,     NULL,    uscript::Compiler::PREC_NONE },       // TOKEN_SUPER
    { NULL,     NULL,    uscript::Compiler::PREC_NONE },       // TOKEN_THIS
    { &uscript::Compiler::Literal,     NULL,    uscript::Compiler::PREC_NONE },       // TOKEN_TRUE
    { NULL,     NULL,    uscript::Compiler::PREC_NONE },       // TOKEN_VAR
    { NULL,     NULL,    uscript::Compiler::PREC_NONE },       // TOKEN_WHILE
    { NULL,     NULL,    uscript::Compiler::PREC_NONE },       // TOKEN_ERROR
    { NULL,     NULL,    uscript::Compiler::PREC_NONE },       // TOKEN_EOF
  };

  return &rules[type];

}

const char *uscript::Compiler::DoIntern(const std::string &str) {
  auto ret = strings.insert(str);
  return ret.first->c_str();
}

void uscript::Compiler::DoRegister(const char *classname) {
  uscript::TClassInfo *ret = tclasslist.Add(classname);
  assert(ret != NULL);
}

