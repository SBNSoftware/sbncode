#include <iostream>

#include "value.h"

void uscript::Value::Print() const {
  switch (val) {
    case uscript::VAL_BOOL: std::cout << (AS_BOOL(*this) ? "true" : "false"); break;
    case uscript::VAL_NIL: std::cout << "nil"; break;
    case uscript::VAL_NUMBER: std::cout << AS_NUMBER(*this); break;
    case uscript::VAL_INTEGER: std::cout << AS_INTEGER(*this); break;
    case uscript::VAL_OBJ_STRING: std::cout << AS_CSTRING(*this); break;
    case uscript::VAL_OBJ_TINSTANCE: std::cout << "TInstance"; break;
    case uscript::VAL_OBJ_TMETHOD: std::cout << "TMethod"; break;
  }
}

bool uscript::Value::operator!() const {
  return IS_NIL(*this) || (IS_INTEGER(*this) && AS_INTEGER(*this) == 0) || (IS_BOOL(*this) && !AS_BOOL(*this));
}

