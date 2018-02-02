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

extern "C" {
  extern ProcessorBase* CreateProcessorObject();
  extern void DestroyProcessorObject(ProcessorBase* proc); 
  extern struct ProcessorBase::export_table exports;
}

/** Load a Processor. */
inline ProcessorBase::export_table* LoadProcessor() {
  return &exports;
  //return (ProcessorBase::export_table *) &exports;
}

int main(int argc, char* argv[]) {
  // Parse command line arguments
  std::vector<char*> config_names;

  int c;
  while ((c=getopt(argc, argv, "c:")) != -1) {
    switch (c) {
      case 'c':
        config_names.push_back(optarg);
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

  if (argc - optind < 1) {
    std::cout << "Usage: " << argv[0] << " [-c [Config]] "
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
  int n_processors = config_names.size() == 0 ? 1 : config_names.size();
  std::vector<Json::Value*> configs(n_processors);
  std::vector<ProcessorBase*> procs(n_processors);

  std::cout << "Configuring... " << std::endl;
  for (size_t i=0; i<n_processors; i++) {
    Json::Value* config = config_names.size() == 0 ? NULL : Main::LoadConfig(config_names[i]);

    procs.push_back(LoadProcessor()->create());
    configs.push_back(config);
  }

  ProcessorBlock block = Main::InitializeBlock(configs, procs);
  std::cout << "Running... " << std::endl;
  block.ProcessFiles(filenames);

  block.DestroyProcessors();
  std::cout << "Done!" << std::endl;

  return 0;
}
