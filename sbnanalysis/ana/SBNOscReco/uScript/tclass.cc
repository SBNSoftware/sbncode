#include <algorithm>

#include "tclass.h"
#include "compile.h"

#include "TClassTable.h"
#include "TClass.h"
#include "TDictionary.h"
#include "TDictAttributeMap.h"
#include "TProtoClass.h"
#include "TSystem.h"
#include "TParameter.h"
#include "TFunction.h"

uscript::TClassInfo *uscript::TClassList::Add(const char *classname) {
  classname = uscript::Compiler::Intern(std::string(classname));
  auto search = classes.find(classname);
  // already added -- ok
  if (search != classes.end()) {
    return &search->second;
  }

  std::map<const char *, uscript::TField> fields;

  DictFuncPtr_t classdict = gClassTable->GetDict(classname);
  if (!classdict) return NULL;

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
      this_field.info = NULL; // no tclass info
      this_field.type = uscript::FIELD_UNSIGNED;
      assert(member->GetUnitSize() == 4);// make sure size is actually 4
    }
    // non-enum basic type
    else if (basic_type != NULL) {
      this_field.info = NULL; // no tclass info
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
      uscript::TClassInfo *field_info = Add(member->GetTypeName());
      if (field_info != NULL) {
        this_field.type = uscript::FIELD_TINSTANCE;
        this_field.info = field_info;
      }
      else {
        can_use = false;
      }
    }
    if (can_use) {
      fields[uscript::Compiler::Intern(std::string(obj->GetName()))] = this_field;
    }
  }
  uscript::TClassInfo tclassinfo;
  tclassinfo.fields = std::move(fields);
  tclassinfo.name = classname;
  return &classes.insert({classname, tclassinfo}).first->second;
}

