#ifndef uscript_value_h
#define uscript_value_h
#include <string>

namespace uscript {

struct Obj;

enum ValueType {
  VAL_BOOL,
  VAL_NIL,
  VAL_NUMBER,
  VAL_INTEGER,
  VAL_OBJ
};

struct Value {
  ValueType val;
  union {
    bool boolean;
    double number;
    int integer;
    Obj *obj;
  } as;
  void Print() const;
  void PrintObj() const;

 bool operator!() const;
};
} // end namespace

#define BOOL_VAL(value)    ((uscript::Value){ uscript::VAL_BOOL, { .boolean = value }})
#define NIL_VAL            ((uscript::Value){ uscript::VAL_NIL,  { .integer = 0 }})
#define NUMBER_VAL(value)  ((uscript::Value){ uscript::VAL_NUMBER, { .number = value }})
#define INTEGER_VAL(value) ((uscript::Value){ uscript::VAL_INTEGER, { .integer = value }})
#define OBJ_VAL(value)     ((uscript::Value){ uscript::VAL_OBJ, { .obj = value }})

#define AS_BOOL(value)    ((value).as.boolean)
#define AS_NUMBER(value)  ((value).as.number)
#define AS_INTEGER(value) ((value).as.integer)
#define AS_OBJ(value)     ((value).as.obj)

#define IS_NIL(value)     ((value).val == uscript::VAL_NIL)
#define IS_BOOL(value)    ((value).val == uscript::VAL_BOOL)
#define IS_NUMBER(value)  ((value).val == uscript::VAL_NUMBER)
#define IS_INTEGER(value) ((value).val == uscript::VAL_INTEGER)
#define IS_OBJ(value)     ((value).val == uscript::VAL_OBJ)

#endif
