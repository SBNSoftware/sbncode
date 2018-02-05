#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <dlfcn.h>
#include <json/json.h>
#include <core/ProcessorBase.hh>
#include <core/ProcessorBlock.hh>

namespace core {

export_table* LoadProcessor(char* libname) {
  std::cout << "Loading processor: " << libname << "... ";

  char* libdir = getenv("SBN_LIB_DIR");
  char libpath[1000];
  snprintf(libpath, 1000, "%s/libsbnanalysis_%s.so", libdir, libname);

  void* handle = dlopen(libpath, RTLD_NOW);
  if (!handle) {
    std::cout << "ERROR" << std::endl;
    std::cerr << "Error loading processor " << libname << ": \n\n"
              << dlerror()
              << "\n\nSBN_LIB_DIR = " << libdir
              << std::endl;
    exit(1);
  }

  export_table* exports = (export_table*) dlsym(handle, "exports");

  if (exports) {
    std::cout << "OK" << std::endl;
  }
  else {
    std::cout << "ERROR" << std::endl;
    std::cerr << "Invalid library " << libpath << ". "
              << "Not a Processor?" << std::endl;
    exit(1);
  }

  return exports;
}


Json::Value* LoadConfig(char* configfile) {
  Json::Value* config = NULL;

  if (configfile) {
    std::cout << "Loading configuration: " << configfile << "... ";
    std::ifstream configstream(configfile, std::ifstream::binary);
    config = new Json::Value;
    Json::Reader reader;
    bool r = reader.parse(configstream, *(config));
    if (!r) {
      std::cerr << "Error parsing configuration file "
                << configfile << std::endl;
      exit(2);
    }
    std::cout << "OK" << std::endl;
  }

  return config;
}

}  // namespace core

