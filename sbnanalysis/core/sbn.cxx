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


/** Load a Processor. */
export_table* LoadProcessor(char* libname) {
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

  export_table* exports = (export_table*) dlsym(handle, "exports");

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


/** Load JSON configuration file if provided. */
Json::Value* LoadConfig(char* configfile) {
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


int main(int argc, char* argv[]) {
  // Parse command line arguments
  std::vector<char*> processors;
  std::map<unsigned, char*> configs;

  int c;
  unsigned procindex = 0;
  while ((c=getopt(argc, argv, "m:c:")) != -1) {
    switch (c) {
      case 'm':
        processors.push_back(optarg);
        procindex++;
        break;
      case 'c':
        configs[procindex-1] = optarg;
        break;
      case '?':
        if (optopt == 'c' || optopt == 'm')
          fprintf(stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint(optopt))
          fprintf(stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
        return 1;
      default:
        abort();
    }
  }

  if (argc - optind < 1) {
    std::cout << "Usage: " << argv[0] << " [-m PROCESSOR [-c CONFIG]] "
              << "INPUTDEF [...]" << std::endl;
    return 0;
  }

  // Process input file definition
  std::vector<std::string> filenames;
  for (int i=optind; i<argc; i++) {
    filenames.push_back(argv[i]);
  }
  assert(!filenames.empty());

  // Setup
  core::ProcessorBlock block;
  std::vector<export_table*> exports(processors.size());
  std::vector<core::ProcessorBase*> procs(processors.size());

  std::cout << "Configuring... " << std::endl;
  for (size_t i=0; i<processors.size(); i++) {
    export_table* exp = LoadProcessor(processors[i]);
    exports[i] = exp;

    Json::Value* config = LoadConfig(configs[i]);

    core::ProcessorBase* proc = exp->create();
    procs[i] = proc;

    block.AddProcessor(proc, config);
  }

  std::cout << "Running... " << std::endl;
  block.ProcessFiles(filenames);

  for (size_t i=0; i<processors.size(); i++) {
    exports[i]->destroy(procs[i]);
  }

  std::cout << "Done!" << std::endl;

  return 0;
}

