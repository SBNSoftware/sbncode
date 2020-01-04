#include <iostream>
#include <string.h>

#include "scanner.h"

uscript::Scanner::Scanner():
  start(NULL),
  current(NULL)
{}

void uscript::Scanner::SetSource(const char *_start) {
  start = _start;
  current = _start;
}

void uscript::Scanner::SkipWhitespace() {
  while (1) {
    char c = *current;
    switch (c) {
      case '/': // comments
        if (current[1] == '/') {
          while (*current != '\n' && !IsAtEnd()) Advance();
        }
        else return;
        break;
      case '\n':
      case ' ':
      case '\r':
      case '\t':
        Advance();
        break;
      default:
        return;
    }
  }
}

static bool isDigit(char c) {
  return c >= '0' && c <= '9';
}

static bool isAlpha(char c) {
  return (c >= 'a' && c <= 'z') ||
         (c >= 'A' && c <= 'Z') ||
          c == '_';
}

uscript::Token uscript::Scanner::ScanToken() {
  SkipWhitespace();

  start = current;

  if (IsAtEnd()) {
    return MakeToken(uscript::TOKEN_EOF);
  }

  char c = Advance();

  if (isDigit(c)) return Number();
  if (isAlpha(c)) return Identifier();

  switch (c) {
    case '(': return MakeToken(uscript::TOKEN_LEFT_PAREN);
    case ')': return MakeToken(uscript::TOKEN_RIGHT_PAREN);
    case '{': return MakeToken(uscript::TOKEN_LEFT_BRACE);
    case '}': return MakeToken(uscript::TOKEN_RIGHT_BRACE);
    case ';': return MakeToken(uscript::TOKEN_SEMICOLON);
    case ',': return MakeToken(uscript::TOKEN_COMMA);
    case '.': return MakeToken(uscript::TOKEN_DOT);
    case '-': return MakeToken(uscript::TOKEN_MINUS);
    case '+': return MakeToken(uscript::TOKEN_PLUS);
    case '/': return MakeToken(uscript::TOKEN_SLASH);
    case '*': return MakeToken(uscript::TOKEN_STAR);
    case '!':
      return MakeToken(Match('=') ? TOKEN_BANG_EQUAL : TOKEN_BANG);
    case '=':
      return MakeToken(Match('=') ? TOKEN_EQUAL_EQUAL : TOKEN_EQUAL);
    case '<':
      return MakeToken(Match('=') ? TOKEN_LESS_EQUAL : TOKEN_LESS);
    case '>':
      return MakeToken(Match('=') ?
                       TOKEN_GREATER_EQUAL : TOKEN_GREATER);
    case '"': return String();
  }

  return ErrorToken("Unexpected character.");
}

bool uscript::Scanner::Match(char c) {
  if (IsAtEnd()) return false;
  if (*current != c) return false;

  current ++;
  return true;
}

char uscript::Scanner::Advance() {
  current++;
  return current[-1];
}

bool uscript::Scanner::IsAtEnd() const {
  return *current == '\0';
}

uscript::Token uscript::Scanner::MakeToken(uscript::TokenType type) const {
  uscript::Token token;
  token.type = type; 
  token.start = start;
  token.length = current - start;
  return token;
}

uscript::Token uscript::Scanner::String() {
  while (*current != '"' && !IsAtEnd()) {
    Advance();
  }

  if (IsAtEnd()) return ErrorToken("Unterminated string.");

  Advance();
  return MakeToken(uscript::TOKEN_STRING);
}

uscript::Token uscript::Scanner::Number() {
  while (isDigit(*current)) Advance();

  if (*current == '.') { // && isDigit(current[1])) {
    // consume the '.'
    Advance();

    while (isDigit(*current)) Advance();
  }

  return MakeToken(uscript::TOKEN_NUMBER);
}

uscript::TokenType uscript::Scanner::IdentifierType() const {
  static const std::vector<std::string> keywords {"and", "or", "if", "else", "for", "while", "return", "true", "false", "fun", "nil", "var", "print"};
  static const std::vector<uscript::TokenType> keyword_types {
    uscript::TOKEN_AND,
    uscript::TOKEN_OR,
    uscript::TOKEN_IF,
    uscript::TOKEN_ELSE,
    uscript::TOKEN_FOR,
    uscript::TOKEN_WHILE,
    uscript::TOKEN_RETURN,
    uscript::TOKEN_TRUE,
    uscript::TOKEN_FALSE,
    uscript::TOKEN_FUN,
    uscript::TOKEN_NIL,
    uscript::TOKEN_VAR,
    uscript::TOKEN_PRINT
  };
  std::string comp = std::string(start, current);
  for (unsigned i = 0; i < keywords.size(); i++) {
    if (strcmp(comp.c_str(), keywords[i].c_str()) == 0) {
      return keyword_types[i];
    }
  }

  return uscript::TOKEN_IDENTIFIER;
}

uscript::Token uscript::Scanner::Identifier() {
  while (isAlpha(*current) || isDigit(*current)) Advance();

  return MakeToken(IdentifierType());
}

uscript::Token uscript::Scanner::ErrorToken(const char *message) const {
  uscript::Token token;
  token.type = uscript::TOKEN_ERROR; 
  token.start = message;
  token.length = (int)strlen(message);
  return token;
}
