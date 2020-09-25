#ifndef uscript_scanner_h
#define uscript_scanner_h

#include "vm.h"
#include "chunk.h"

namespace uscript {
enum TokenType {
  // Single-character tokens.
  TOKEN_LEFT_PAREN, TOKEN_RIGHT_PAREN,
  TOKEN_LEFT_BRACE, TOKEN_RIGHT_BRACE,
  TOKEN_LEFT_BRACKET, TOKEN_RIGHT_BRACKET,
  TOKEN_COMMA, TOKEN_DOT, TOKEN_MINUS, TOKEN_PLUS,
  TOKEN_SEMICOLON, TOKEN_SLASH, TOKEN_STAR,

  // One or two character tokens.
  TOKEN_BANG, TOKEN_BANG_EQUAL,
  TOKEN_EQUAL, TOKEN_EQUAL_EQUAL,
  TOKEN_GREATER, TOKEN_GREATER_EQUAL,
  TOKEN_LESS, TOKEN_LESS_EQUAL,

  // Literals.
  TOKEN_IDENTIFIER, TOKEN_STRING, TOKEN_NUMBER,

  // Keywords.
  TOKEN_AND, TOKEN_CLASS, TOKEN_ELSE, TOKEN_FALSE,
  TOKEN_FOR, TOKEN_FUN, TOKEN_IF, TOKEN_NIL, TOKEN_OR,
  TOKEN_PRINT, TOKEN_RETURN, TOKEN_SUPER, TOKEN_THIS,
  TOKEN_TRUE, TOKEN_VAR, TOKEN_WHILE,
  
  // more keywords
  TOKEN_LENGTH, TOKEN_FIELDS,

  TOKEN_ERROR,
  TOKEN_EOF
};

struct Token {
  TokenType type;
  const char *start;
  int length;
};

class Scanner {
  const char *start;
  const char *current;
public:
  Scanner();
  void SetSource(const char *_start);

  Token ScanToken();
  Token MakeToken(TokenType type) const;
  Token ErrorToken(const char *message) const;
  bool IsAtEnd() const;
  char Advance();
  bool Match(char c);
  void SkipWhitespace();
  Token String();
  Token Number();
  Token Identifier();
  TokenType IdentifierType() const;
};


} // end namespace
#endif
