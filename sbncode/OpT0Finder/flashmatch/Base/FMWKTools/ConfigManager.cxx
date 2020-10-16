#ifndef __FLASHMATCHBASE_CONFIGMANAGER_CXX__
#define __FLASHMATCHBASE_CONFIGMANAGER_CXX__

#include "ConfigManager.h"
#include "PSetUtils.h"
namespace flashmatch {

  ConfigManager* ConfigManager::_me = nullptr;

  void ConfigManager::AddConfigFile(const std::string cfg_file)
  {
    if(_cfg_files.find(cfg_file)!=_cfg_files.end()) {
      std::cerr << "Duplicate addition of config fiel: " << cfg_file << std::endl;
      throw std::exception();
    }

    _cfg.add_pset(CreatePSetFromFile(cfg_file));
  }

  void ConfigManager::AddConfigString(const std::string cfg_str)
  {
    PSet p;
    p.add(cfg_str);
    _cfg.add_pset(p);
  }
  
  const PSet& ConfigManager::GetConfig(const std::string cfg)
  {
    return _cfg.get_pset(cfg);
  }

}

#endif
