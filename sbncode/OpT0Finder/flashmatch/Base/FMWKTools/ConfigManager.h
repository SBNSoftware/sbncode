/**
 * \file ConfigManager.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class flashmatch::ConfigManager
 *
 * @author drinkingkazu
 */

/** \addtogroup core_Base

    @{*/
#ifndef __FLASHMATCHBASE_CONFIGMANAGER_H__
#define __FLASHMATCHBASE_CONFIGMANAGER_H__

#include <iostream>
#include "PSet.h"
#include <set>

namespace flashmatch {
  /**
     \class ConfigManager
     \brief Utility class to register a set of configurations
     Provides also a shared instance through which registered configurations can be shared beyond a single owner.\n
     Using flashmatch::PSet, the uniqueness of configuration parameters is guaranteed (no worry to "overwrite")\n
  */
  class ConfigManager {
    
  public:
    
    /// Default constructor
    ConfigManager() {}
         
    /// Default destructor
    ~ConfigManager(){}
    /// Shared static reference getter
    static const ConfigManager& get() 
    {
      if(!_me) _me = new ConfigManager;
      return *_me;
    }
    /// Adder of configuration from a file
    void AddConfigFile(const std::string cfg_file);
    /// Adder of configuration from parsed string
    void AddConfigString(const std::string cfg_str);
    /// Configuration retrieval method
    const PSet& GetConfig(const std::string cfg);

  private:

    static ConfigManager* _me;
    std::set<std::string> _cfg_files;
    PSet _cfg;
    
  };
}
#endif
/** @} */ // end of doxygen group 

