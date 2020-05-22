#ifndef uscript_value_h
#define uscript_value_h
#include <string>

#include "tclass.h"

namespace uscript {

enum ValueType {
  VAL_BOOL,
  VAL_NIL,
  VAL_NUMBER,
  VAL_INTEGER,
  VAL_OBJ_STRING,
  VAL_OBJ_TINSTANCE,
};

struct ObjString {
  const char *string;
};

struct TClassInfo;

struct ObjTInstance {
  uint8_t *loc;
  TData data; 
};

struct Value {
  ValueType val;
  union {
    bool boolean;
    double number;
    int integer;
    ObjString string;
    ObjTInstance tinst;
  } as;
  void Print() const;
  void PrintObj() const;

 bool operator!() const;
};
} // end namespace

#define BOOL_VAL(value)      ((uscript::Value){ uscript::VAL_BOOL, { .boolean = value }})
#define NIL_VAL              ((uscript::Value){ uscript::VAL_NIL,  { .integer = 0 }})
#define NUMBER_VAL(value)    ((uscript::Value){ uscript::VAL_NUMBER, { .number = value }})
#define INTEGER_VAL(value)   ((uscript::Value){ uscript::VAL_INTEGER, { .integer = value }})
#define STRING_VAL(value)    ((uscript::Value){ uscript::VAL_OBJ_STRING, { .string = (uscript::ObjString){ value } }})
#define TINSTANCE_VAL(value) ((uscript::Value){ uscript::VAL_OBJ_TINSTANCE, { .tinst = value }})

#define AS_BOOL(value)    ((value).as.boolean)
#define AS_NUMBER(value)  ((value).as.number)
#define AS_INTEGER(value) ((value).as.integer)
#define AS_STRING(value)  ((value).as.string)
#define AS_CSTRING(value) ((value).as.string.string)
#define AS_TINSTANCE(value) ((value).as.tinst)

#define IS_NIL(value)       ((value).val == uscript::VAL_NIL)
#define IS_BOOL(value)      ((value).val == uscript::VAL_BOOL)
#define IS_NUMBER(value)    ((value).val == uscript::VAL_NUMBER)
#define IS_INTEGER(value)   ((value).val == uscript::VAL_INTEGER)
#define IS_STRING(value)    ((value).val == uscript::VAL_OBJ_STRING)
#define IS_TINSTANCE(value) ((value).val == uscript::VAL_OBJ_TINSTANCE)

#endif
