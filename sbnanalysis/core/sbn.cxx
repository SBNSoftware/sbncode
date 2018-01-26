#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <dlfcn.h>
#include <json/json.h>
#include <core/ProcessorBase.hh>

int main(int argc, char* argv[]) {
  // Command line arguments
  char* configfile = NULL;
  int c;

  while ((c=getopt(argc, argv, "c:")) != -1) {
    switch (c) {
      case 'c':
        configfile = optarg;
        break;
      case '?':
        if (optopt == 'c')
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

  if (argc - optind < 2) {
    std::cout << "Usage: " << argv[0] << " [-c CONFIG] "
              << "PROCESSOR INPUTDEF [...]" << std::endl;
    return 0;
  }

  // Load JSON configuration file if provided
  Json::Value* config = NULL;
  if (configfile) {
    std::cout << "Configuration: " << configfile << std::endl;

    std::ifstream configstream(configfile, std::ifstream::binary);
    config = new Json::Value;
    Json::Reader reader;
    bool r = reader.parse(configstream, *config);
    assert(r);
  }

  // Process input file definition
  std::vector<std::string> filenames;
  for (int i=optind+1; i<argc; i++) {
    filenames.push_back(argv[i]);
  }
  assert(!filenames.empty());

  // Load the requested Processor plugin
  char* libname = argv[optind];
  std::cout << "Processor: " << libname << std::endl;

  char* libdir = getenv("SBN_LIB_DIR");
  char libpath[1000];
  snprintf(libpath, 1000, "%s/libsbnanalysis_%s.so", libdir, libname);

  std::cout << "Loading processor " << libname << "..." << std::endl;

  void* handle = dlopen("lib/libsbnanalysis_ExampleAnalysis_ExampleSelection.so", RTLD_LAZY);
  if (!handle) {
    std::cerr << "Processor " << libname << " not found in SBN_LIB_DIR " << libdir << std::endl;
    return 1;
  }

  // Set up pointers to Processor create/delete functions
  core::ProcessorBase* (*create)();
  void (*destroy)(core::ProcessorBase*);

  create = (core::ProcessorBase* (*)()) dlsym(handle, "CreateObject");
  destroy = (void (*)(core::ProcessorBase*)) dlsym(handle, "DestroyObject");

  // Run the Processor
  core::ProcessorBase* proc = (core::ProcessorBase*) create();
  std::cout << "Running " << libname << "..." << std::endl;
  proc->ProcessFiles(filenames, config);
  destroy(proc);

  std::cout << "Done!" << std::endl;

  return 0;
}

