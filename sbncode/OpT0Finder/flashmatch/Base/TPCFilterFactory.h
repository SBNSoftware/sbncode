/**
 * \file TPCFilterFactory.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class TPCFilterFactory
 *
 * @author kazuhiro
 */

/** \addtogroup Base

    @{*/
#ifndef TPCFILTERFACTORY_H
#define TPCFILTERFACTORY_H

#include <iostream>
#include <map>
#include "BaseTPCFilter.h"
namespace flashmatch {

  /**
     \class TPCFilterFactoryBase
     \brief Abstract base class for factory (to be implemented per flash)
  */
  class TPCFilterFactoryBase {
  public:
    /// Default ctor
    TPCFilterFactoryBase(){}
    /// Default dtor (virtual)
    virtual ~TPCFilterFactoryBase(){}
    /// Abstract constructor method
    virtual BaseTPCFilter* create(const std::string instance_name) = 0;
  };

  /**
     \class TPCFilterFactory
     \brief Factory class for instantiating flash algorithm instance
  */
  class TPCFilterFactory {
  private:
    /// Default ctor, shouldn't be used
    TPCFilterFactory() {}
  public:
    /// Default dtor
    ~TPCFilterFactory() {_factory_map.clear();}
    /// Static sharable instance getter
    static TPCFilterFactory& get()
    { if(!_me) _me = new TPCFilterFactory; return *_me; }
    /// Factory registration method (should be called by global factory instance in algorithm header)
    void add_factory(const std::string name, flashmatch::TPCFilterFactoryBase* factory)
    { _factory_map[name] = factory; }
    /// Factory creation method (should be called by clients, possibly you!)
    BaseTPCFilter* create(const std::string name, const std::string instance_name) {
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
    std::map<std::string,flashmatch::TPCFilterFactoryBase*> _factory_map;
    /// Static self
    static TPCFilterFactory* _me;
  };
}
#endif
/** @} */ // end of doxygen group 

