/**
 * \file ChargeAnalytical.h
 *
 * \ingroup Algorithms
 * 
 * \brief Class def header for a class ChargeAnalytical
 *
 * @author matthias
 */

/** \addtogroup Algorithms

    @{*/
#ifndef CHARGEANALYTICAL_H
#define CHARGEANALYTICAL_H

#ifndef USING_LARSOFT
#define USING_LARSOFT 1
#endif

#if USING_LARSOFT == 0
#include "flashmatch/Base/BaseFlashFilter.h"
#include "flashmatch/Base/FlashHypothesisFactory.h"
#include "flashmatch/Base/FMWKInterface.h"
#include "flashmatch/Base/OpT0FinderException.h"
#else
#include "sbncode/OpT0Finder/flashmatch/Base/BaseFlashFilter.h"
#include "sbncode/OpT0Finder/flashmatch/Base/FlashHypothesisFactory.h"
#include "sbncode/OpT0Finder/flashmatch/Base/FMWKInterface.h"
#include "sbncode/OpT0Finder/flashmatch/Base/OpT0FinderException.h"
#endif

#include <iostream>
#include <math.h>

namespace flashmatch {

  /**
     \class ChargeAnalytical
     User defined class ChargeAnalytical ... these comments are used to generate
     doxygen documentation!
  */
  class ChargeAnalytical : public BaseFlashHypothesis {
    
  public:
    
    /// Default constructor
    ChargeAnalytical(const std::string name="ChargeAnalytical");
    
    /// Default destructor
    ~ChargeAnalytical() { }
    
    void FillEstimate(const QCluster_t&, Flash_t&) const;

  protected:

    void _Configure_(const Config_t &pset);

    double _global_qe;         ///< Global QE
    std::vector<double> _qe_v; ///< PMT-wise relative QE
    
  };

  /**
     \class flashmatch::ChargeAnalyticalFactory
  */
  class ChargeAnalyticalFactory : public FlashHypothesisFactoryBase {
  public:
    /// ctor
    ChargeAnalyticalFactory() { FlashHypothesisFactory::get().add_factory("ChargeAnalytical",this); }
    /// dtor
    ~ChargeAnalyticalFactory() {}
    /// creation method
    BaseFlashHypothesis* create(const std::string instance_name) { return new ChargeAnalytical(instance_name); }
  };
  
}

#endif
/** @} */ // end of doxygen group 

