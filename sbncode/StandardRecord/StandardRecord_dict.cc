// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME StandardRecord_dict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "classes.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *cafcLcLSRHeader_Dictionary();
   static void cafcLcLSRHeader_TClassManip(TClass*);
   static void *new_cafcLcLSRHeader(void *p = 0);
   static void *newArray_cafcLcLSRHeader(Long_t size, void *p);
   static void delete_cafcLcLSRHeader(void *p);
   static void deleteArray_cafcLcLSRHeader(void *p);
   static void destruct_cafcLcLSRHeader(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::caf::SRHeader*)
   {
      ::caf::SRHeader *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::caf::SRHeader));
      static ::ROOT::TGenericClassInfo 
         instance("caf::SRHeader", 28, "SRHeader.h", 15,
                  typeid(::caf::SRHeader), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &cafcLcLSRHeader_Dictionary, isa_proxy, 12,
                  sizeof(::caf::SRHeader) );
      instance.SetNew(&new_cafcLcLSRHeader);
      instance.SetNewArray(&newArray_cafcLcLSRHeader);
      instance.SetDelete(&delete_cafcLcLSRHeader);
      instance.SetDeleteArray(&deleteArray_cafcLcLSRHeader);
      instance.SetDestructor(&destruct_cafcLcLSRHeader);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::caf::SRHeader*)
   {
      return GenerateInitInstanceLocal((::caf::SRHeader*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::caf::SRHeader*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *cafcLcLSRHeader_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::caf::SRHeader*)0x0)->GetClass();
      cafcLcLSRHeader_TClassManip(theClass);
   return theClass;
   }

   static void cafcLcLSRHeader_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *cafcLcLSRSlice_Dictionary();
   static void cafcLcLSRSlice_TClassManip(TClass*);
   static void *new_cafcLcLSRSlice(void *p = 0);
   static void *newArray_cafcLcLSRSlice(Long_t size, void *p);
   static void delete_cafcLcLSRSlice(void *p);
   static void deleteArray_cafcLcLSRSlice(void *p);
   static void destruct_cafcLcLSRSlice(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::caf::SRSlice*)
   {
      ::caf::SRSlice *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::caf::SRSlice));
      static ::ROOT::TGenericClassInfo 
         instance("caf::SRSlice", 23, "SRSlice.h", 11,
                  typeid(::caf::SRSlice), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &cafcLcLSRSlice_Dictionary, isa_proxy, 12,
                  sizeof(::caf::SRSlice) );
      instance.SetNew(&new_cafcLcLSRSlice);
      instance.SetNewArray(&newArray_cafcLcLSRSlice);
      instance.SetDelete(&delete_cafcLcLSRSlice);
      instance.SetDeleteArray(&deleteArray_cafcLcLSRSlice);
      instance.SetDestructor(&destruct_cafcLcLSRSlice);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::caf::SRSlice*)
   {
      return GenerateInitInstanceLocal((::caf::SRSlice*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::caf::SRSlice*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *cafcLcLSRSlice_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::caf::SRSlice*)0x0)->GetClass();
      cafcLcLSRSlice_TClassManip(theClass);
   return theClass;
   }

   static void cafcLcLSRSlice_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *cafcLcLSRSpill_Dictionary();
   static void cafcLcLSRSpill_TClassManip(TClass*);
   static void *new_cafcLcLSRSpill(void *p = 0);
   static void *newArray_cafcLcLSRSpill(Long_t size, void *p);
   static void delete_cafcLcLSRSpill(void *p);
   static void deleteArray_cafcLcLSRSpill(void *p);
   static void destruct_cafcLcLSRSpill(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::caf::SRSpill*)
   {
      ::caf::SRSpill *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::caf::SRSpill));
      static ::ROOT::TGenericClassInfo 
         instance("caf::SRSpill", 36, "SRSpill.h", 15,
                  typeid(::caf::SRSpill), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &cafcLcLSRSpill_Dictionary, isa_proxy, 12,
                  sizeof(::caf::SRSpill) );
      instance.SetNew(&new_cafcLcLSRSpill);
      instance.SetNewArray(&newArray_cafcLcLSRSpill);
      instance.SetDelete(&delete_cafcLcLSRSpill);
      instance.SetDeleteArray(&deleteArray_cafcLcLSRSpill);
      instance.SetDestructor(&destruct_cafcLcLSRSpill);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::caf::SRSpill*)
   {
      return GenerateInitInstanceLocal((::caf::SRSpill*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::caf::SRSpill*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *cafcLcLSRSpill_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::caf::SRSpill*)0x0)->GetClass();
      cafcLcLSRSpill_TClassManip(theClass);
   return theClass;
   }

   static void cafcLcLSRSpill_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *cafcLcLStandardRecord_Dictionary();
   static void cafcLcLStandardRecord_TClassManip(TClass*);
   static void *new_cafcLcLStandardRecord(void *p = 0);
   static void *newArray_cafcLcLStandardRecord(Long_t size, void *p);
   static void delete_cafcLcLStandardRecord(void *p);
   static void deleteArray_cafcLcLStandardRecord(void *p);
   static void destruct_cafcLcLStandardRecord(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::caf::StandardRecord*)
   {
      ::caf::StandardRecord *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::caf::StandardRecord));
      static ::ROOT::TGenericClassInfo 
         instance("caf::StandardRecord", 29, "StandardRecord.h", 24,
                  typeid(::caf::StandardRecord), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &cafcLcLStandardRecord_Dictionary, isa_proxy, 12,
                  sizeof(::caf::StandardRecord) );
      instance.SetNew(&new_cafcLcLStandardRecord);
      instance.SetNewArray(&newArray_cafcLcLStandardRecord);
      instance.SetDelete(&delete_cafcLcLStandardRecord);
      instance.SetDeleteArray(&deleteArray_cafcLcLStandardRecord);
      instance.SetDestructor(&destruct_cafcLcLStandardRecord);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::caf::StandardRecord*)
   {
      return GenerateInitInstanceLocal((::caf::StandardRecord*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::caf::StandardRecord*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *cafcLcLStandardRecord_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::caf::StandardRecord*)0x0)->GetClass();
      cafcLcLStandardRecord_TClassManip(theClass);
   return theClass;
   }

   static void cafcLcLStandardRecord_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_cafcLcLSRHeader(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::caf::SRHeader : new ::caf::SRHeader;
   }
   static void *newArray_cafcLcLSRHeader(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::caf::SRHeader[nElements] : new ::caf::SRHeader[nElements];
   }
   // Wrapper around operator delete
   static void delete_cafcLcLSRHeader(void *p) {
      delete ((::caf::SRHeader*)p);
   }
   static void deleteArray_cafcLcLSRHeader(void *p) {
      delete [] ((::caf::SRHeader*)p);
   }
   static void destruct_cafcLcLSRHeader(void *p) {
      typedef ::caf::SRHeader current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::caf::SRHeader

namespace ROOT {
   // Wrappers around operator new
   static void *new_cafcLcLSRSlice(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::caf::SRSlice : new ::caf::SRSlice;
   }
   static void *newArray_cafcLcLSRSlice(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::caf::SRSlice[nElements] : new ::caf::SRSlice[nElements];
   }
   // Wrapper around operator delete
   static void delete_cafcLcLSRSlice(void *p) {
      delete ((::caf::SRSlice*)p);
   }
   static void deleteArray_cafcLcLSRSlice(void *p) {
      delete [] ((::caf::SRSlice*)p);
   }
   static void destruct_cafcLcLSRSlice(void *p) {
      typedef ::caf::SRSlice current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::caf::SRSlice

namespace ROOT {
   // Wrappers around operator new
   static void *new_cafcLcLSRSpill(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::caf::SRSpill : new ::caf::SRSpill;
   }
   static void *newArray_cafcLcLSRSpill(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::caf::SRSpill[nElements] : new ::caf::SRSpill[nElements];
   }
   // Wrapper around operator delete
   static void delete_cafcLcLSRSpill(void *p) {
      delete ((::caf::SRSpill*)p);
   }
   static void deleteArray_cafcLcLSRSpill(void *p) {
      delete [] ((::caf::SRSpill*)p);
   }
   static void destruct_cafcLcLSRSpill(void *p) {
      typedef ::caf::SRSpill current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::caf::SRSpill

namespace ROOT {
   // Wrappers around operator new
   static void *new_cafcLcLStandardRecord(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::caf::StandardRecord : new ::caf::StandardRecord;
   }
   static void *newArray_cafcLcLStandardRecord(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::caf::StandardRecord[nElements] : new ::caf::StandardRecord[nElements];
   }
   // Wrapper around operator delete
   static void delete_cafcLcLStandardRecord(void *p) {
      delete ((::caf::StandardRecord*)p);
   }
   static void deleteArray_cafcLcLStandardRecord(void *p) {
      delete [] ((::caf::StandardRecord*)p);
   }
   static void destruct_cafcLcLStandardRecord(void *p) {
      typedef ::caf::StandardRecord current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::caf::StandardRecord

namespace {
  void TriggerDictionaryInitialization_libStandardRecord_dict_Impl() {
    static const char* headers[] = {
0    };
    static const char* includePaths[] = {
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "libStandardRecord_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace caf{class __attribute__((annotate("$clingAutoload$SRHeader.h")))  __attribute__((annotate("$clingAutoload$StandardRecord.h")))  SRHeader;}
namespace caf{class __attribute__((annotate("$clingAutoload$SRSlice.h")))  __attribute__((annotate("$clingAutoload$StandardRecord.h")))  SRSlice;}
// namespace caf{class __attribute__((annotate("$clingAutoload$SRSpill.h")))  __attribute__((annotate("$clingAutoload$StandardRecord.h")))  SRSpill;}
namespace caf{class __attribute__((annotate("$clingAutoload$StandardRecord.h")))  StandardRecord;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "libStandardRecord_dict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif
#ifndef USE_ROOT
  #define USE_ROOT 1
#endif
#ifndef _POSIX_SOURCE
  #define _POSIX_SOURCE 1
#endif
#ifndef _SVID_SOURCE
  #define _SVID_SOURCE 1
#endif
#ifndef _BSD_SOURCE
  #define _BSD_SOURCE 1
#endif
#ifndef _POSIX_C_SOURCE
  #define _POSIX_C_SOURCE 2
#endif
#ifndef _XOPEN_SOURCE
  #define _XOPEN_SOURCE 1
#endif
#ifndef UNIX
  #define UNIX 1
#endif
#ifndef LINUX
  #define LINUX 1
#endif
#ifndef __UNIX__
  #define __UNIX__ 1
#endif
#ifndef __LINUX__
  #define __LINUX__ 1
#endif
#ifndef DATAREP_LITTLE_IEEE
  #define DATAREP_LITTLE_IEEE 1
#endif
#ifndef DATAREP_LITTLE_ENDIAN
  #define DATAREP_LITTLE_ENDIAN 1
#endif
#ifndef DEFECT_NO_IOSTREAM_NAMESPACES
  #define DEFECT_NO_IOSTREAM_NAMESPACES 1
#endif
#ifndef DEFECT_NO_JZEXT
  #define DEFECT_NO_JZEXT 1
#endif
#ifndef DEFECT_NO_INTHEX
  #define DEFECT_NO_INTHEX 1
#endif
#ifndef DEFECT_NO_INTHOLLERITH
  #define DEFECT_NO_INTHOLLERITH 1
#endif
#ifndef DEFECT_NO_READONLY
  #define DEFECT_NO_READONLY 1
#endif
#ifndef DEFECT_NO_DIRECT_FIXED
  #define DEFECT_NO_DIRECT_FIXED 1
#endif
#ifndef DEFECT_NO_STRUCTURE
  #define DEFECT_NO_STRUCTURE 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "StandardRecord.h"

// #include "SRSpillTruthBranch.h" // not part of the main record

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"caf::SRHeader", payloadCode, "@",
"caf::SRSlice", payloadCode, "@",
"caf::SRSpill", payloadCode, "@",
"caf::StandardRecord", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libStandardRecord_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libStandardRecord_dict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libStandardRecord_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libStandardRecord_dict() {
  TriggerDictionaryInitialization_libStandardRecord_dict_Impl();
}
