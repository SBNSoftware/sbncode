/**
 * \file MatchSelectionFactory.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class MatchSelectionFactory
 *
 * @author kazuhiro
 */

/** \addtogroup Base

    @{*/
#ifndef MatchSelectionFACTORY_H
#define MatchSelectionFACTORY_H

#include <iostream>
#include <map>
#include "BaseMatchSelection.h"
namespace flashmatch {

  /**
     \class MatchSelectionFactoryBase
     \brief Abstract base class for factory (to be implemented per flash)
  */
  class MatchSelectionFactoryBase {
  public:
    /// Default ctor
    MatchSelectionFactoryBase(){}
    /// Default dtor (virtual)
    virtual ~MatchSelectionFactoryBase(){}
    /// Abstract constructor method
    virtual BaseMatchSelection* create(const std::string instance_name) = 0;
  };

  /**
     \class MatchSelectionFactory
     \brief Factory class for instantiating flash algorithm instance
  */
  class MatchSelectionFactory {
  private:
    /// Default ctor, shouldn't be used
    MatchSelectionFactory() {}
  public:
    /// Default dtor
    ~MatchSelectionFactory() {_factory_map.clear();}
    /// Static sharable instance getter
    static MatchSelectionFactory& get()
    { if(!_me) _me = new MatchSelectionFactory; return *_me; }
    /// Factory registration method (should be called by global factory instance in algorithm header)
    void add_factory(const std::string name, flashmatch::MatchSelectionFactoryBase* factory)
    { _factory_map[name] = factory; }
    /// Factory creation method (should be called by clients, possibly you!)
    BaseMatchSelection* create(const std::string name, const std::string instance_name) {
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
    std::map<std::string,flashmatch::MatchSelectionFactoryBase*> _factory_map;
    /// Static self
    static MatchSelectionFactory* _me;
  };
}
#endif
/** @} */ // end of doxygen group 

