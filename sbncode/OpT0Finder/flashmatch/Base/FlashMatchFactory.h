/**
 * \file FlashMatchFactory.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class FlashMatchFactory
 *
 * @author kazuhiro
 */

/** \addtogroup Base

    @{*/
#ifndef FLASHMATCHFACTORY_H
#define FLASHMATCHFACTORY_H

#include <iostream>
#include <map>
#include "BaseFlashMatch.h"
namespace flashmatch {

  /**
     \class FlashMatchFactoryBase
     \brief Abstract base class for factory (to be implemented per flash)
  */
  class FlashMatchFactoryBase {
  public:
    /// Default ctor
    FlashMatchFactoryBase(){}
    /// Default dtor (virtual)
    virtual ~FlashMatchFactoryBase(){}
    /// Abstract constructor method
    virtual BaseFlashMatch* create(const std::string instance_name) = 0;
  };

  /**
     \class FlashMatchFactory
     \brief Factory class for instantiating flash algorithm instance
  */
  class FlashMatchFactory {
  private:
    /// Default ctor, shouldn't be used
    FlashMatchFactory() {}
  public:
    /// Default dtor
    ~FlashMatchFactory() {_factory_map.clear();}
    /// Static sharable instance getter
    static FlashMatchFactory& get()
    { if(!_me) _me = new FlashMatchFactory; return *_me; }
    /// Factory registration method (should be called by global factory instance in algorithm header)
    void add_factory(const std::string name, flashmatch::FlashMatchFactoryBase* factory)
    { _factory_map[name] = factory; }
    /// Factory creation method (should be called by clients, possibly you!)
    BaseFlashMatch* create(const std::string name, const std::string instance_name) {
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
    std::map<std::string,flashmatch::FlashMatchFactoryBase*> _factory_map;
    /// Static self
    static FlashMatchFactory* _me;
  };
}
#endif
/** @} */ // end of doxygen group 

