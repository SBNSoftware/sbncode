#include <iostream>
#include <string>

#include <sys/stat.h>

#include "TError.h"
#include "TFile.h"
#include "TTree.h"

int main(int argc, char** argv)
{
  gErrorIgnoreLevel = 100000;

  if(argc != 2){
    std::cerr << "Usage: must supply one filename as an argument" << std::endl;
    exit(1);
  }

  const std::string filePath = argv[1];

  struct stat buf;
  if(stat(filePath.c_str(), &buf) != 0){
    std::cerr << "ERROR: File does not exist: " << filePath << std::endl;
    exit(1);
  }

  std::unique_ptr<TFile> f(TFile::Open(filePath.c_str(), "READ"));

  if(!f->IsOpen()){
    std::cerr << "ERROR: Unable to open " << filePath
              << " as a TFile, is this a proper ROOT file?" << std::endl;
    exit(1);
  }

  TDirectory* metadata = f->GetDirectory("metadata");
  if(!metadata){
    std::cerr << "ERROR: Unable to access metadata in " << filePath
              << " is this a proper CAF with metadata?" << std::endl;
    exit(1);
  }

  TTree* tr = (TTree*)metadata->Get("metatree");
  if(!tr){
    std::cerr << "ERROR: Unable to access metadata tree in " << filePath
              << " is this a proper CAF with metadata?" << std::endl;
    exit(1);
  }

  std::cout << "{\n";

  std::string key, value;
  std::string* pkey = &key;
  std::string* pvalue = &value;
  tr->SetBranchAddress("key", &pkey);
  tr->SetBranchAddress("value", &pvalue);

  bool first = true;
  for(int i = 0; i < tr->GetEntries(); ++i){
    if(!first) std::cout << "," << std::endl;
    first = false;

    tr->GetEntry(i);

    std::cout << "  \"" << key << "\": ";
    //      if(value.empty()) value = "\"none\"";
    // The convention is that values are already suitably escaped inside the
    // CAF file.
    std::cout << value;
  }
  std::cout << "\n}" << std::endl;

  return 0;
}
