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


extern "C" {
  extern core::ProcessorBase* CreateProcessorObject();
  extern void DestroyProcessorObject(core::ProcessorBase* proc); 
  extern struct processor_export_table processor_exports;
}

/** Load a Processor. */
export_table* LoadProcessor() {
  return (export_table *) &processor_exports;
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
  std::map<unsigned, char*> configs;

  int c;
  unsigned procindex = 0;
  while ((c=getopt(argc, argv, "c:")) != -1) {
    switch (c) {
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
  std::vector<export_table*> exports(configs.size());
  std::vector<core::ProcessorBase*> procs(configs.size());

  std::cout << "Configuring... " << std::endl;
  for (size_t i=0; i<procs.size(); i++) {
    export_table* exp = LoadProcessor();
    exports[i] = exp;

    Json::Value* config = LoadConfig(configs[i]);

    core::ProcessorBase* proc = exp->create();
    procs[i] = proc;

    block.AddProcessor(proc, config);
  }

  std::cout << "Running... " << std::endl;
  block.ProcessFiles(filenames);

  for (size_t i=0; i<procs.size(); i++) {
    exports[i]->destroy(procs[i]);
  }

  std::cout << "Done!" << std::endl;

  return 0;
}

