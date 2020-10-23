#ifndef __PSET_UTILS_CXX__
#define __PSET_UTILS_CXX__

#include "PSetUtils.h"
#include <sstream>
#include <fstream>
namespace flashmatch {

  std::string ConfigFile2String(std::string fname)
  {
    std::ifstream filestrm(fname.c_str());
    std::string   contents;
    std::string   line;

    while(std::getline(filestrm, line)) {

      if(line.empty()) continue;
      
      std::stringstream linestrm(line);
      std::string       valid_line;
      
      std::getline(linestrm, valid_line, '#');
      
      if(valid_line.empty()) continue;
      
      contents += " " + valid_line;
    }
    filestrm.close();    
    return contents;
  }

  PSet CreatePSetFromFile(std::string fname,std::string cfg_name)
  {
    PSet res(cfg_name,ConfigFile2String(fname));
    return res;
  }

}

#endif
