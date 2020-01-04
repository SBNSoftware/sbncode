#include <algorithm>

#include "tclass.h"

#include "TClassTable.h"
#include "TClass.h"
#include "TDictionary.h"
#include "TDictAttributeMap.h"
#include "TProtoClass.h"
#include "TSystem.h"
#include "TParameter.h"
#include "TFunction.h"

int uscript::TClassList::Add(const char *classname) {
  int index = std::distance(classnames.begin(), std::find(classnames.begin(), classnames.end(), std::string(classname)));
  // already added -- ok
  if (index < classnames.size()) {
    return index;
  }

  std::map<std::string, uscript::TField> fields;

  DictFuncPtr_t classdict = gClassTable->GetDict(classname);
  if (!classdict) return -1;

  TClass *tclass = classdict();

  TList *members = tclass->GetListOfDataMembers();
  TIterator *m_iterator = members->MakeIterator();
  TObject *obj;
  while ((obj = m_iterator->Next()) != NULL) {
    TDataMember *member = (TDataMember *)obj;

    uscript::TField this_field;
    this_field.offset = member->GetOffset();
    bool can_use = true;
    TDataType *basic_type = member->GetDataType();

    // treat enum as unsigned
    if (member->IsEnum()) {
      this_field.tclassIndex = -1; // no tclass index      
      this_field.type = uscript::FIELD_UNSIGNED;
      assert(member->GetUnitSize() == 4);// make sure size is actually 4
    }
    // non-enum basic type
    else if (basic_type != NULL) {
      this_field.tclassIndex = -1; // no tclass index
      switch (basic_type->GetType()) {
        case kBool_t:
          this_field.type = uscript::FIELD_BOOL;
          break;
        case kFloat_t:
          this_field.type = uscript::FIELD_FLOAT;
          break;
        case kDouble_t:
          this_field.type = uscript::FIELD_DOUBLE;
          break;
        case kInt_t:
          this_field.type = uscript::FIELD_INT;
          break;
        case kUInt_t:
          this_field.type = uscript::FIELD_UNSIGNED;
          break;
        default: // don't allow any other field type
          can_use = false;
          break;
      }
    }
    // complex-type (TClass)
    else {
      // must have a class dictionary
      int index = Add(member->GetTypeName());
      if (index >= 0) {
        this_field.type = uscript::FIELD_TINSTANCE;
        this_field.tclassIndex = index;
      }
      else {
        can_use = false;
      }
    }
    if (can_use) {
      fields[std::string(obj->GetName())] = this_field;
    }
  }
  uscript::TClassInfo tclassinfo;
  tclassinfo.fields = std::move(fields);
  tclassinfo.name = classname;
  classes.push_back(std::move(tclassinfo));
  return classes.size() - 1;
}

