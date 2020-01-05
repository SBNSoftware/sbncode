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
#include "TClassEdit.h"

uscript::TClassInfo *uscript::TClassList::Add(const char *classname) {
  classname = uscript::Compiler::Intern(std::string(classname));
  auto search = classes.find(classname);
  // already added -- ok
  if (search != classes.end()) {
    return &search->second;
  }

  DictFuncPtr_t classdict = gClassTable->GetDict(classname);

  // no classdict -- bad
  if (!classdict) return NULL;

  std::map<const char *, uscript::TField> fields;

  TClass *tclass = classdict();

  // handle vector case
  if (tclass->GetCollectionType() == ROOT::kSTLvector) {
    TClassInfo tclassinfo;
    tclassinfo.is_vec = true;
    tclassinfo.name = classname; 
    tclassinfo.size = sizeof(std::vector<void*>); // all vectors in layout should have the same size

    TClassEdit::TSplitType split(classname);
    // get the dictionary for the sub-class
    classdict = gClassTable->GetDict(split.fElements[1].c_str());

    // vec data
    uscript::TData data;
    data.len = -1;

    if (classdict == NULL) {
      // try to see if it is a basic type
      if (split.fElements[1] == "bool" || split.fElements[1] == "Bool_t") {
        data.type = uscript::FIELD_BOOL; 
      }
      else if (split.fElements[1] == "float" || split.fElements[1] == "Float_t") {
        data.type = uscript::FIELD_FLOAT;
      }
      else if (split.fElements[1] == "double" || split.fElements[1] == "Double_t") {
        data.type = uscript::FIELD_DOUBLE;
      }
      else if (split.fElements[1] == "int" || split.fElements[1] == "Int_t") {
        data.type = uscript::FIELD_INT;
      }
      else if (split.fElements[1] == "unsigned" || split.fElements[1] == "UInt_t") {
        data.type = uscript::FIELD_UNSIGNED;
      }
      else return NULL;
    }
    // a dictionary exists -- try to add it
    else {
      TClassInfo *sub_info = Add(split.fElements[1].c_str());
      data.info = sub_info;
      data.type = uscript::FIELD_TINSTANCE;
    }

    tclassinfo.vec_data = data;
    return &classes.insert({classname, tclassinfo}).first->second;
  }

  TList *members = tclass->GetListOfDataMembers();
  TIterator *m_iterator = members->MakeIterator();
  TObject *obj;
  while ((obj = m_iterator->Next()) != NULL) {
    TDataMember *member = (TDataMember *)obj;

    uscript::TField this_field;
    this_field.offset = member->GetOffset();
    bool can_use = true;
    TDataType *basic_type = member->GetDataType();

    // get size
    this_field.data.len = -1;
    if (member->GetArrayDim() > 0) {
      this_field.data.len = 1;
      for (unsigned i = 0; i < member->GetArrayDim(); i++) {
        this_field.data.len = this_field.data.len * member->GetMaxIndex(i);
      }
    }

    // treat enum as unsigned with size 4
    if (member->IsEnum()) {
      this_field.data.info = NULL; // no tclass info
      this_field.data.type = uscript::FIELD_ENUM;
      assert(member->GetUnitSize() == 4);// make sure size is actually 4
    }
    // non-enum basic type
    else if (basic_type != NULL) {
      this_field.data.info = NULL; // no tclass info
      switch (basic_type->GetType()) {
        case kBool_t:
          this_field.data.type = uscript::FIELD_BOOL;
          break;
        case kFloat_t:
          this_field.data.type = uscript::FIELD_FLOAT;
          break;
        case kDouble_t:
          this_field.data.type = uscript::FIELD_DOUBLE;
          break;
        case kInt_t:
          this_field.data.type = uscript::FIELD_INT;
          break;
        case kUInt_t:
          this_field.data.type = uscript::FIELD_UNSIGNED;
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
        this_field.data.type = uscript::FIELD_TINSTANCE;
        this_field.data.info = field_info;
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
  tclassinfo.size = tclass->GetClassSize();
  tclassinfo.is_vec = false;
  return &classes.insert({classname, tclassinfo}).first->second;
}

int uscript::TData::Size() const {
  switch (type) {
    case uscript::FIELD_BOOL: return sizeof(bool);
    case uscript::FIELD_INT:  return sizeof(int);
    case uscript::FIELD_UNSIGNED: return sizeof(unsigned);
    case uscript::FIELD_FLOAT: return sizeof(float);
    case uscript::FIELD_DOUBLE: return sizeof(double);
    case uscript::FIELD_ENUM: return 4;
    case uscript::FIELD_TINSTANCE: return info->size;
  }
  // unreachable
  return -1;
}

int uscript::TData::Length(uint8_t *loc) const {
  if (info && info->is_vec) {
    // FIXME: this is probably implementation defined -- needs to be fixed
    uint8_t **vec = (uint8_t**)loc;
    uint8_t *start = vec[0];
    uint8_t *end = vec[1];
    return (end - start) / Size();
  }
  else {
    return len;
  }

}

