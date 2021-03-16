#include <iostream>
#include <string>

#include <sys/stat.h>

#include "TError.h"
#include "TFile.h"
#include "TTree.h"

std::string escapeQuotes(std::string str)
{
  size_t pos = 0;
  while(true){
    pos = str.find("\"", pos);
    if(pos == std::string::npos) return str;
    str.replace(pos, 1, "\\\"");
    pos += 2; // Where to start looking next time
  }
}

int main(int argc, char** argv)
{
  gErrorIgnoreLevel = 100000;

  if(argc < 2){
    std::cerr << "Usage: must supply at least one filename as an argument"
              << std::endl;
    exit(1);
  }

  std::cout << "{" << std::endl;

  for(int i = 1; i < argc; ++i){
    const std::string filePath = argv[i];

    struct stat buf;
    if(stat(filePath.c_str(), &buf) != 0){
      std::cerr << "ERROR: File does not exist: " << filePath << std::endl;
      continue;
    }
    
    const char* base = basename(filePath.c_str());

    std::unique_ptr<TFile> f(TFile::Open(filePath.c_str(), "READ"));
    
    if(!f->IsOpen()){
      std::cerr << "ERROR: Unable to open " << filePath
                << " as a TFile, is this a proper ROOT file?" << std::endl;
      continue;
    }

    std::cout << "  \"" << base << "\": {" << std::endl;

    TDirectory* metadata = f->GetDirectory("metadata");
    if(!metadata){
      std::cerr << "ERROR: Unable to access metadata in " << filePath
                << " is this a proper CAF with metadata?" << std::endl;
      continue;
    }

    TTree* tr = (TTree*)metadata->Get("metatree");
    if(!tr){
      std::cerr << "ERROR: Unable to access metadata tree in " << filePath
                << " is this a proper CAF with metadata?" << std::endl;
      continue;
    }

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

      std::cout << "    \"" << key << "\": ";
      //      if(value.empty()) value = "none";
      std::cout << "\"" << escapeQuotes(value) << "\"";
    }
    std::cout << "\n  }";
    if(i < argc-1) std::cout << ",";
    std::cout << std::endl;
  } // end for i

  std::cout << "}" << std::endl;

  return 0;
}
