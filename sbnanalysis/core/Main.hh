#ifndef __sbnanalysis_core_MAIN__
#define __sbnanalysis_core_MAIN__

#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <vector>
#include <dlfcn.h>
#include <json/json.h>
#include <core/ProcessorBase.hh>
#include <core/ProcessorBlock.hh>


namespace core {
class Main {
  public:

  struct export_table {
    core::ProcessorBase* (*create)(void);
    void (*destroy)(core::ProcessorBase*);
  };


  struct ProcessorDef {
    char* name;
    struct export_table* exports;
    Json::Value* config;
    core::ProcessorBase* proc;
  };

  /** Load JSON configuration file if provided. */
  static Json::Value* LoadConfig(char* configfile);

  static core::ProcessorBlock InitializeBlock(std::vector<Json::Value*> &export_tables, std::vector<core::ProcessorBase*> &processors); 

  static export_table* LoadProcessor(char* libname);
};
}

#endif
