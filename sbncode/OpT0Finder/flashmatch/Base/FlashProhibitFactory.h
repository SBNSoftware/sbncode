/**
 * \file FlashProhibitFactory.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class FlashProhibitFactory
 *
 * @author kazuhiro
 */

/** \addtogroup Base

    @{*/
#ifndef FLASHPROHIBITFACTORY_H
#define FLASHPROHIBITFACTORY_H

#include <iostream>
#include <map>
#include "BaseProhibitAlgo.h"
namespace flashmatch {

  /**
     \class FlashProhibitFactoryBase
     \brief Abstract base class for factory (to be implemented per flash)
  */
  class FlashProhibitFactoryBase {
  public:
    /// Default ctor
    FlashProhibitFactoryBase(){}
    /// Default dtor (virtual)
    virtual ~FlashProhibitFactoryBase(){}
    /// Abstract constructor method
    virtual BaseProhibitAlgo* create(const std::string instance_name) = 0;
  };

  /**
     \class FlashProhibitFactory
     \brief Factory class for instantiating flash algorithm instance
  */
  class FlashProhibitFactory {
  private:
    /// Default ctor, shouldn't be used
    FlashProhibitFactory() {}
  public:
    /// Default dtor
    ~FlashProhibitFactory() {_factory_map.clear();}
    /// Static sharable instance getter
    static FlashProhibitFactory& get()
    { if(!_me) _me = new FlashProhibitFactory; return *_me; }
    /// Factory registration method (should be called by global factory instance in algorithm header)
    void add_factory(const std::string name, flashmatch::FlashProhibitFactoryBase* factory)
    { _factory_map[name] = factory; }
    /// Factory creation method (should be called by clients, possibly you!)
    BaseProhibitAlgo* create(const std::string name, const std::string instance_name) {
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
    std::map<std::string,flashmatch::FlashProhibitFactoryBase*> _factory_map;
    /// Static self
    static FlashProhibitFactory* _me;
  };
}
#endif
/** @} */ // end of doxygen group 

