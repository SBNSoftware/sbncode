#ifndef uscript_tclass_h
#define uscript_tclass_h

#include <vector>
#include <string>
#include <map>
#include <unordered_map>

namespace uscript {

enum TFieldType {
  FIELD_BOOL,
  FIELD_INT,
  FIELD_UNSIGNED,
  FIELD_FLOAT,
  FIELD_DOUBLE,
  FIELD_ENUM,
  FIELD_TINSTANCE
};

struct TClassInfo;

struct TData {
  TClassInfo *info;
  TFieldType type;
  int len;
  int Size() const;
  int Length(uint8_t *loc) const;
};

struct TField {
  int offset;
  TData data;
};

struct TClassInfo {
  std::map<const char *, TField> fields;
  std::string name;
  bool is_vec;
  TData vec_data;
  int size;
};

class TClassList {
public:
  TClassInfo *Add(const char *classname);
  std::unordered_map<const char *, TClassInfo> classes;

};

} // end namespace
#endif
