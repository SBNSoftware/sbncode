#ifndef uscript_tclass_h
#define uscript_tclass_h

#include <vector>
#include <string>
#include <map>
#include <set>

namespace uscript {

enum TFieldType {
  FIELD_BOOL,
  FIELD_INT,
  FIELD_UNSIGNED,
  FIELD_FLOAT,
  FIELD_DOUBLE,
  FIELD_TINSTANCE
};

struct TField {
  int offset;
  int tclassIndex;
  TFieldType type;
};

struct TClassInfo {
  std::map<std::string, TField> fields;  
  std::string name;
};

class TClassList {
public:
  int Add(const char *classname);
  std::vector<TClassInfo> classes;
  std::set<std::string> classnames;

};

} // end namespace
#endif
