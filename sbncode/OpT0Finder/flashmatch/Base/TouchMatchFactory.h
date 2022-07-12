/**
 * \file TouchMatchFactory.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class TouchMatchFactory
 *
 * @author kazuhiro
 */

/** \addtogroup Base

    @{*/
#ifndef TouchMatchFACTORY_H
#define TouchMatchFACTORY_H

#include <iostream>
#include <map>
#include "BaseTouchMatch.h"
namespace flashmatch {

  /**
     \class TouchMatchFactoryBase
     \brief Abstract base class for factory (to be implemented per flash)
  */
  class TouchMatchFactoryBase {
  public:
    /// Default ctor
    TouchMatchFactoryBase(){}
    /// Default dtor (virtual)
    virtual ~TouchMatchFactoryBase(){}
    /// Abstract constructor method
    virtual BaseTouchMatch* create(const std::string instance_name) = 0;
  };

  /**
     \class TouchMatchFactory
     \brief Factory class for instantiating flash algorithm instance
  */
  class TouchMatchFactory {
  private:
    /// Default ctor, shouldn't be used
    TouchMatchFactory() {}
  public:
    /// Default dtor
    ~TouchMatchFactory() {_factory_map.clear();}
    /// Static sharable instance getter
    static TouchMatchFactory& get()
    { if(!_me) _me = new TouchMatchFactory; return *_me; }
    /// Factory registration method (should be called by global factory instance in algorithm header)
    void add_factory(const std::string name, flashmatch::TouchMatchFactoryBase* factory)
    { _factory_map[name] = factory; }
    /// Factory creation method (should be called by clients, possibly you!)
    BaseTouchMatch* create(const std::string name, const std::string instance_name) {
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
    std::map<std::string,flashmatch::TouchMatchFactoryBase*> _factory_map;
    /// Static self
    static TouchMatchFactory* _me;
  };
}
#endif
/** @} */ // end of doxygen group 

