#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <vector>
#include <dlfcn.h>
#include <fhiclcpp/ParameterSet.h>
#include <core/ProcessorBase.hh>
#include <core/ProcessorBlock.hh>
#include <core/Loader.hh>

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
	break;
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
  std::string filedef = argv[optind];
  std::string list_suffix = ".list";
  std::vector<std::string> filenames;

  if (std::equal(list_suffix.rbegin(), list_suffix.rend(), filedef.rbegin())) {
    // File list
    std::ifstream infile(filedef);
    std::string filename;
    while (infile >> filename) {
      filenames.push_back(filename);
    }
  }
  else {
    // Files listed on command line
    for (int i=optind; i<argc; i++) {
      filenames.push_back(argv[i]);
    }
  }

  assert(!filenames.empty());

  // Setup
  std::vector<core::ProcessorBase*> procs(processors.size());
  std::vector<fhicl::ParameterSet*> configs(processors.size());

  std::cout << "Configuring... " << std::endl;
  for (size_t i=0; i<procs.size(); i++) {
    core::export_table* exp = core::LoadProcessor(processors[i]);
    procs[i] = exp->create();
    configs[i] = core::LoadConfig(config_names[i]);
  }

  core::ProcessorBlock block;
  for (int i=0; i<procs.size(); i++) {
    block.AddProcessor(procs[i], configs[i]);
  }

  std::cout << "Running... " << std::endl;
  block.ProcessFiles(filenames);

  block.DeleteProcessors();
  std::cout << "Done!" << std::endl;

  return 0;
}

