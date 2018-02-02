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
#include <core/Main.hh>

using namespace core;

/** Load JSON configuration file if provided. */
Json::Value* Main::LoadConfig(char* configfile) {
  Json::Value* config = NULL;
  
  if (configfile) {
    std::cout << "Loading configuration: " << configfile << "... ";
    std::ifstream configstream(configfile, std::ifstream::binary);
    config = new Json::Value;
    Json::Reader reader;
    bool r = reader.parse(configstream, *(config));
    if (!r) {
      std::cerr << "Error parsing configuration file " << configfile << std::endl;
      exit(2);
    }
    std::cout << "OK" << std::endl;
  }

  return config;
}

ProcessorBlock Main::InitializeBlock(std::vector<Json::Value*> &configs, std::vector<ProcessorBase*> &processors) {
  ProcessorBlock block;
  assert(configs.size() == processors.size());
  int size = configs.size();
  for (int i =0; i < size; i++) {
    block.AddProcessor(processors[i], configs[i]);
  }
  return block;
} 

/** Load a Processor. */
ProcessorBase::export_table* Main::LoadProcessor(char* libname) {
  std::cout << "Loading processor: " << libname << "... ";

  char* libdir = getenv("SBN_LIB_DIR");
  char libpath[1000];
  snprintf(libpath, 1000, "%s/libsbnanalysis_%s.so", libdir, libname);

  void* handle = dlopen(libpath, RTLD_NOW);
  if (!handle) {
    std::cout << "ERROR" << std::endl;
    std::cerr << "Error loading processor " << libpath << ": \n\n"
              << dlerror()
              << "\n\nSBN_LIB_DIR = " << libdir
              << std::endl;
    exit(1);
  }

  ProcessorBase::export_table* exports = (ProcessorBase::export_table*) dlsym(handle, "exports");

  if (exports) {
    std::cout << "OK" << std::endl;
  }
  else {
    std::cout << "ERROR" << std::endl;
    std::cerr << "Invalid library " << libname << ". Not a Processor?" << std::endl;
    exit(1);
  }

  return exports;
}


