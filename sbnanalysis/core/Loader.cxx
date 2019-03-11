#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <dlfcn.h>
#include <cetlib/filepath_maker.h>
#include <fhiclcpp/make_ParameterSet.h>
#include <fhiclcpp/ParameterSet.h>
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

export_table_postprocess* LoadPostProcessor(char* libname) {
  std::cout << "Loading post-processor: " << libname << "... ";

  char* libdir = getenv("SBN_LIB_DIR");
  char libpath[1000];
  snprintf(libpath, 1000, "%s/libsbnanalysis_%s.so", libdir, libname);

  void* handle = dlopen(libpath, RTLD_NOW);
  if (!handle) {
    std::cout << "ERROR" << std::endl;
    std::cerr << "Error loading post-processor " << libname << ": \n\n"
              << dlerror()
              << "\n\nSBN_LIB_DIR = " << libdir
              << std::endl;
    exit(1);
  }

  export_table_postprocess* exports = (export_table_postprocess*) dlsym(handle, "exports");

  if (exports) {
    std::cout << "OK" << std::endl;
  }
  else {
    std::cout << "ERROR" << std::endl;
    std::cerr << "Invalid library " << libpath << ". "
              << "Not a PostProcessor?" << std::endl;
    exit(1);
  }

  return exports;
}


fhicl::ParameterSet* LoadConfig(char* configfile) {
  fhicl::ParameterSet* config = NULL;

  if (configfile) {
    std::cout << "Loading configuration: " << configfile << "... ";
    config = new fhicl::ParameterSet;
    cet::filepath_lookup_nonabsolute maker("FHICL_FILE_PATH");
    fhicl::make_ParameterSet(configfile, maker, *config);
    std::cout << "OK" << std::endl;
  }

  return config;
}

}  // namespace core

