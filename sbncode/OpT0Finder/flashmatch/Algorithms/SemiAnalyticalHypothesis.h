/**
 * \file SemiAnalyticalHypothesis.h
 *
 * \ingroup Algorithms
 *
 * \brief Class def header for a class SemiAnalyticalHypothesis
 *
 * @author mdeltutto
 */

/** \addtogroup Algorithms

    @{*/

#ifndef SEMIANALYTICALHYPOTHESIS_H
#define SEMIANALYTICALHYPOTHESIS_H

#include <iostream>

#ifndef USING_LARSOFT
#define USING_LARSOFT 1
#endif

#include "PhotonLibHypothesis.h"

#if USING_LARSOFT == 0
#include "flashmatch/Base/FlashHypothesisFactory.h"
#include "flashmatch/Base/FMWKInterface.h"
#include "flashmatch/Base/OpT0FinderException.h"
#else
#include "sbncode/OpT0Finder/flashmatch/Base/FlashHypothesisFactory.h"
#include "sbncode/OpT0Finder/flashmatch/Base/FMWKInterface.h"
#include "sbncode/OpT0Finder/flashmatch/Base/OpT0FinderException.h"
#include "larsim/LegacyLArG4/OpFastScintillation.hh"
#endif

namespace flashmatch {
  /**
     \class SemiAnalyticalHypothesis
     User custom analysis class made by SHELL_USER_NAME
   */
  class SemiAnalyticalHypothesis : public PhotonLibHypothesis {

  public:

    /// Default constructor
    SemiAnalyticalHypothesis(const std::string name="SemiAnalyticalHypothesis");

    /// Default destructor
    virtual ~SemiAnalyticalHypothesis(){}

    void BuildHypothesis(const QCluster_t&, Flash_t&) const;

  protected:

    #if USING_LARSOFT == 1
    larg4::OpFastScintillation* _opfast_scintillation; ///< For SBND semi-analytical
    #endif
  };

  /**
     \class flashmatch::SemiAnalyticalHypothesisFactory
  */
  class SemiAnalyticalHypothesisFactory : public FlashHypothesisFactoryBase {
  public:
    /// ctor
    SemiAnalyticalHypothesisFactory() { FlashHypothesisFactory::get().add_factory("SemiAnalyticalHypothesis",this); }
    /// dtor
    ~SemiAnalyticalHypothesisFactory() {}
    /// creation method
    BaseFlashHypothesis* create(const std::string instance_name) { return new SemiAnalyticalHypothesis(instance_name); }
  };
}
#endif

/** @} */ // end of doxygen group
