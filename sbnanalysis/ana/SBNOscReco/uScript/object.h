#ifndef uscript_object_h
#define uscript_object_h

#include <string>

#include "value.h"

namespace uscript {

enum ObjType {
  OBJ_STRING,
  OBJ_TINSTANCE,
  OBJ_TMETHOD
};

struct Obj {
  ObjType type;
};

struct ObjString : public Obj {
  std::string string;
  explicit ObjString(std::string s);
};

struct ObjTInstance: public Obj {
  uint8_t *data;
  int tclassIndex;
};

struct ObjTMethod: public Obj {

};

} // end namespace

#define OBJ_TYPE(value) (AS_OBJ(value)->type)

static inline bool isObjType(uscript::Value value, uscript::ObjType type) {
  return IS_OBJ(value) && AS_OBJ(value)->type == type;
}

#define IS_STRING(value) isObjType(value, uscript::OBJ_STRING)
#define AS_STRING(value) (((uscript::ObjString *)AS_OBJ(value))->string)
#define AS_CSTRING(value) (((uscript::ObjString *)AS_OBJ(value))->string.c_str())
#define AS_TMETHOD(value) ((uscript::ObjTMethod *)AS_OBJ(value))
#define AS_TINSTANCE(value) ((uscript::ObjTInstance *)AS_OBJ(value))

#define IS_TINSTANCE(value) isObjType(value, uscript::OBJ_TINSTANCE)

#endif
