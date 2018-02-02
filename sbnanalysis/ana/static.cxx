#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <dlfcn.h>
#include <json/json.h>
#include <core/ProcessorBase.hh>
#include CMAKE_PROCESSOR_INCLUDE

int main(int argc, char* argv[]) {
  // Parse command line arguments
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

  /*
  if (argc - optind < 2) {
    std::cout << "Usage: " << argv[0] << " [-c CONFIG] "
              << "PROCESSOR INPUTDEF [...]" << std::endl;
    return 0;
  }*/

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

  void (*destroy)(core::ProcessorBase*);
  destroy = (void (*)(core::ProcessorBase*)) DestroyObject;

  // Run the Processor
  core::ProcessorBase* proc = (core::ProcessorBase*) CreateObject();
  proc->ProcessFiles(filenames, config);
  destroy(proc);

  std::cout << "Done!" << std::endl;

  return 0;
}

