/**
 * An example analysis.
 */

#include <string>
#include <vector>
#include "ExampleTools.h"
#include "ExampleSelection.h"

using namespace ana::ExampleAnalysis;

int main(int argc, char* argv[]) {
  std::vector<std::string> filenames = { argv[1] };

  // Use a library function (from ExampleTools)
  hello();

  // Set up an analysis Processor
  ExampleSelection proc;

  proc.Initialize();

  proc.ProcessFile(filenames);

  proc.Finalize();

  return 0;
}

