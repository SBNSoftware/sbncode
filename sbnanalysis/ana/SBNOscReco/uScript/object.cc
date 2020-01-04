#include "object.h"

uscript::ObjString::ObjString(std::string s):
  string(s)
{
  type = uscript::OBJ_STRING;
}

