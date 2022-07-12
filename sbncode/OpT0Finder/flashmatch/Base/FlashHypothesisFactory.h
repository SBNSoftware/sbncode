/**
 * \file FlashHypothesisFactory.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class FlashHypothesisFactory
 *
 * @author kazuhiro
 */

/** \addtogroup Base

    @{*/
#ifndef FLASHHYPOTHESISFACTORY_H
#define FLASHHYPOTHESISFACTORY_H

#include <iostream>
#include <map>
#include "BaseFlashHypothesis.h"
namespace flashmatch {

  /**
     \class FlashHypothesisFactoryBase
     \brief Abstract base class for factory (to be implemented per flash)
  */
  class FlashHypothesisFactoryBase {
  public:
    /// Default ctor
    FlashHypothesisFactoryBase(){}
    /// Default dtor (virtual)
    virtual ~FlashHypothesisFactoryBase(){}
    /// Abstract constructor method
    virtual BaseFlashHypothesis* create(const std::string instance_name) = 0;
  };

  /**
     \class FlashHypothesisFactory
     \brief Factory class for instantiating flash algorithm instance
  */
  class FlashHypothesisFactory {
  private:
    /// Default ctor, shouldn't be used
    FlashHypothesisFactory() {}
  public:
    /// Default dtor
    ~FlashHypothesisFactory() {_factory_map.clear();}
    /// Static sharable instance getter
    static FlashHypothesisFactory& get()
    { if(!_me) _me = new FlashHypothesisFactory; return *_me; }
    /// Factory registration method (should be called by global factory instance in algorithm header)
    void add_factory(const std::string name, flashmatch::FlashHypothesisFactoryBase* factory)
    { _factory_map[name] = factory; }
    /// Factory creation method (should be called by clients, possibly you!)
    BaseFlashHypothesis* create(const std::string name, const std::string instance_name) {
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
    std::map<std::string,flashmatch::FlashHypothesisFactoryBase*> _factory_map;
    /// Static self
    static FlashHypothesisFactory* _me;
  };
}
#endif
/** @} */ // end of doxygen group 

