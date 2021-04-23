/**
 * \file BaseAlgorithm.h
 *
 * \ingroup Base
 * 
 * \brief Class def header for a class BaseAlgorithm
 *
 * @author kazuhiro
 */

/** \addtogroup Base

    @{*/
#ifndef OPT0FINDER_BASEALGORITHM_H
#define OPT0FINDER_BASEALGORITHM_H

#include "OpT0FinderTypes.h"
#include "FMWKInterface.h"
#include "LoggerFeature.h"
#include <vector>
namespace flashmatch {

  class FlashMatchManager;
  
  /**
     \class BaseAlgorithm
  */
  class BaseAlgorithm : public LoggerFeature {

    friend class FlashMatchManager;
    
  public:
    
    /// Default constructor
    BaseAlgorithm(const Algorithm_t type, const std::string name);
    
    /// Default destructor
    ~BaseAlgorithm(){}

    /// Function to accept configuration
    void Configure(const Config_t &pset);

    /// Algorithm type
    Algorithm_t AlgorithmType() const;

    /// Algorithm name
    const std::string& AlgorithmName() const;

  protected:

    virtual void _Configure_(const Config_t &pset) = 0;

  private:
    
    Algorithm_t _type; ///< Algorithm type
    std::string _name; ///< Algorithm name
  };
}
#endif
/** @} */ // end of doxygen group 

