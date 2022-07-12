/**
 * \file CustomAlgoFactory.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class CustomAlgoFactory
 *
 * @author kazuhiro
 */

/** \addtogroup Base

    @{*/
#ifndef CUSTOMFILTERFACTORY_H
#define CUSTOMFILTERFACTORY_H

#include <iostream>
#include <map>
#include "BaseAlgorithm.h"
namespace flashmatch {

  /**
     \class CustomAlgoFactoryBase
     \brief Abstract base class for factory (to be implemented per flash)
  */
  class CustomAlgoFactoryBase {
  public:
    /// Default ctor
    CustomAlgoFactoryBase(){}
    /// Default dtor (virtual)
    virtual ~CustomAlgoFactoryBase(){}
    /// Abstract constructor method
    virtual BaseAlgorithm* create(const std::string instance_name) = 0;
  };

  /**
     \class CustomAlgoFactory
     \brief Factory class for instantiating flash algorithm instance
  */
  class CustomAlgoFactory {
  private:
    /// Default ctor, shouldn't be used
    CustomAlgoFactory() {}
  public:
    /// Default dtor
    ~CustomAlgoFactory() {_factory_map.clear();}
    /// Static sharable instance getter
    static CustomAlgoFactory& get()
    { if(!_me) _me = new CustomAlgoFactory; return *_me; }
    /// Factory registration method (should be called by global factory instance in algorithm header)
    void add_factory(const std::string name, flashmatch::CustomAlgoFactoryBase* factory)
    { _factory_map[name] = factory; }
    /// Factory creation method (should be called by clients, possibly you!)
    BaseAlgorithm* create(const std::string name, const std::string instance_name) {
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
    std::map<std::string,flashmatch::CustomAlgoFactoryBase*> _factory_map;
    /// Static self
    static CustomAlgoFactory* _me;
  };
}
#endif
/** @} */ // end of doxygen group 

