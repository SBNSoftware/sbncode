/**
 * \file FlashFilterFactory.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class FlashFilterFactory
 *
 * @author kazuhiro
 */

/** \addtogroup Base

    @{*/
#ifndef FLASHFILTERFACTORY_H
#define FLASHFILTERFACTORY_H

#include <iostream>
#include <map>
#include "BaseFlashFilter.h"
namespace flashmatch {

  /**
     \class FlashFilterFactoryBase
     \brief Abstract base class for factory (to be implemented per flash)
  */
  class FlashFilterFactoryBase {
  public:
    /// Default ctor
    FlashFilterFactoryBase(){}
    /// Default dtor (virtual)
    virtual ~FlashFilterFactoryBase(){}
    /// Abstract constructor method
    virtual BaseFlashFilter* create(const std::string instance_name) = 0;
  };

  /**
     \class FlashFilterFactory
     \brief Factory class for instantiating flash algorithm instance
  */
  class FlashFilterFactory {
  private:
    /// Default ctor, shouldn't be used
    FlashFilterFactory() {}
  public:
    /// Default dtor
    ~FlashFilterFactory() {_factory_map.clear();}
    /// Static sharable instance getter
    static FlashFilterFactory& get()
    { if(!_me) _me = new FlashFilterFactory; return *_me; }
    /// Factory registration method (should be called by global factory instance in algorithm header)
    void add_factory(const std::string name, flashmatch::FlashFilterFactoryBase* factory)
    { _factory_map[name] = factory; }
    /// Factory creation method (should be called by clients, possibly you!)
    BaseFlashFilter* create(const std::string name, const std::string instance_name) {
      auto iter = _factory_map.find(name);
      if(iter == _factory_map.end() || !((*iter).second)) {
	std::cerr << "Found no registered class " << name << std::endl;
	return nullptr;
      }
      auto ptr = (*iter).second->create(instance_name);
      return ptr;
    }

  private:
    /// Static factory container
    std::map<std::string,flashmatch::FlashFilterFactoryBase*> _factory_map;
    /// Static self
    static FlashFilterFactory* _me;
  };
}
#endif
/** @} */ // end of doxygen group 

