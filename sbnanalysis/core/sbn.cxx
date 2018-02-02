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

int main(int argc, char* argv[]) {
  // Parse command line arguments
  std::vector<char*> processors;
  std::map<unsigned, char*> config_names;

  int c;
  unsigned procindex = 0;
  while ((c=getopt(argc, argv, "m:c:")) != -1) {
    switch (c) {
      case 'm':
        processors.push_back(optarg);
        procindex++;
        break;
      case 'c':
        config_names[procindex-1] = optarg;
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
  std::vector<ProcessorBase*> procs(processors.size());
  std::vector<Json::Value*> configs(processors.size());

  std::cout << "Configuring... " << std::endl;
  for (size_t i=0; i<processors.size(); i++) {
    Main::export_table* exp = Main::LoadProcessor(processors[i]);

    Json::Value* config = Main::LoadConfig(config_names[i]);
    configs[i] = config;

    ProcessorBase* proc = exp->create();
    procs[i] = proc;
  }

  ProcessorBlock block = Main::InitializeBlock(configs, procs);

  std::cout << "Running... " << std::endl;
  block.ProcessFiles(filenames);

  block.DestroyProcessors();
  std::cout << "Done!" << std::endl;

  return 0;
}

