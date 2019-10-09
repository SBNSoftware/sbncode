#include <fstream>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <vector>
#include <dlfcn.h>
#include <core/PostProcessorBase.hh>
#include <core/Loader.hh>

int main(int argc, char* argv[]) {
  // Parse command line arguments
  char *post_processor = NULL;
  char *config_name = NULL; // optional
  std::string output_fname = "";

  int c;
  unsigned procindex = 0;
  while ((c=getopt(argc, argv, "m:c:o:")) != -1) {
    switch (c) {
      case 'm':
        if (post_processor != NULL) {
          fprintf(stderr, "Only specify one post-processor.\n");
          return 1;
        }
        post_processor = optarg;
        break;
      case 'c':
        if (config_name != NULL) {
          fprintf(stderr, "Only specify one post-processor configuration.\n");
          return 1;
        }
        config_name = optarg;
        break;
      case 'o':
        output_fname = optarg; 
        break;
      case '?':
        if (optopt == 'c' || optopt == 'm' || optopt == 'o')
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
    std::cout << "Usage: " << argv[0] << " -m POSTPROCESSOR [-c CONFIG] "
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
  assert(post_processor != NULL);

  // Setup
  std::cout << "Configuring... " << std::endl;
  core::export_table_postprocess *exp = core::LoadPostProcessor(post_processor);
  core::PostProcessorBase *proc = exp->create();
  proc->Initialize(config_name, output_fname);

  // Run
  std::cout << "Running... " << std::endl;
  proc->Run(filenames);
  std::cout << "Done!" << std::endl;

  return 0;
}

