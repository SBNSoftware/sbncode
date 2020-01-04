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
  FIELD_TINSTANCE
};

struct TClassInfo;

struct TField {
  int offset;
  TClassInfo *info;
  TFieldType type;
};

struct TClassInfo {
  std::map<const char *, TField> fields;
  std::string name;
};

class TClassList {
public:
  TClassInfo *Add(const char *classname);
  std::unordered_map<const char *, TClassInfo> classes;

};

} // end namespace
#endif
