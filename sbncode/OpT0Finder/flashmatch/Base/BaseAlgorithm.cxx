#ifndef OPT0FINDER_BASEALGORITHM_CXX
#define OPT0FINDER_BASEALGORITHM_CXX

#include "BaseAlgorithm.h"
#include "OpT0FinderException.h"
namespace flashmatch {

  BaseAlgorithm::BaseAlgorithm(const Algorithm_t type,const std::string name)
    : LoggerFeature(name)
    , _type(type)
    , _name(name)
  {}

  Algorithm_t BaseAlgorithm::AlgorithmType() const
  { return _type; }

  void BaseAlgorithm::Configure(const Config_t &pset)
  {
    this->set_verbosity((msg::Level_t)(pset.get<unsigned int>("Verbosity",(unsigned int)(msg::kNORMAL))));
    this->_Configure_(pset);
  }

  const std::string& BaseAlgorithm::AlgorithmName() const
  { return _name; }

}

#endif
