#include <iostream>

#include "value.h"
#include "object.h"

void uscript::Value::Print() const {
  switch (val) {
    case uscript::VAL_BOOL: std::cout << (AS_BOOL(*this) ? "true" : "false"); break;
    case uscript::VAL_NIL: std::cout << "nil"; break;
    case uscript::VAL_NUMBER: std::cout << AS_NUMBER(*this) << " num"; break;
    case uscript::VAL_INTEGER: std::cout << AS_INTEGER(*this) << " int"; break;
    case uscript::VAL_OBJ: PrintObj(); break;
  }
}

void uscript::Value::PrintObj() const {
  switch (OBJ_TYPE(*this)) {
    case uscript::OBJ_STRING:
      std::cout << AS_STRING(*this);
      break;
    case uscript::OBJ_TINSTANCE:
      std::cout << "TInstance";
      break;
    default:
      break;
  }
}

bool uscript::Value::operator!() const {
return IS_NIL(*this) || (IS_INTEGER(*this) && AS_INTEGER(*this) == 0) || (IS_BOOL(*this) && !AS_BOOL(*this));
}

