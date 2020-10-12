/**
 * \file PhotonLibHypothesis.h
 *
 * \ingroup Algorithms
 *
 * \brief Class def header for a class PhotonLibHypothesis
 *
 * @author yuntse
 */

/** \addtogroup Algorithms

    @{*/

#ifndef PHOTONLIBHYPOTHESIS_H
#define PHOTONLIBHYPOTHESIS_H

#ifndef USING_LARSOFT
#define USING_LARSOFT 1
#endif

#if USING_LARSOFT == 0
#include "flashmatch/Base/OpT0FinderException.h"
#include "flashmatch/Base/FMWKInterface.h"
#include "flashmatch/Base/BaseFlashFilter.h"
#include "flashmatch/Base/FlashHypothesisFactory.h"
#else
#include "sbncode/OpT0Finder/flashmatch/Base/OpT0FinderException.h"
#include "sbncode/OpT0Finder/flashmatch/Base/FMWKInterface.h"
#include "sbncode/OpT0Finder/flashmatch/Base/BaseFlashFilter.h"
#include "sbncode/OpT0Finder/flashmatch/Base/FlashHypothesisFactory.h"
#include "larsim/LegacyLArG4/OpFastScintillation.hh"
#endif

#include <iostream>
#include <cassert>
#include <chrono>

namespace flashmatch {
  /**
     \class PhotonLibHypothesis
     User custom analysis class made by SHELL_USER_NAME
   */
  class PhotonLibHypothesis : public BaseFlashHypothesis {

  public:

    /// Default constructor
    PhotonLibHypothesis(const std::string name="PhotonLibHypothesis");

    /// Default destructor
    virtual ~PhotonLibHypothesis(){}

    void FillEstimate(const QCluster_t&, Flash_t&) const;

  private:
    /// Fills the estimate using the semi analytical approach (SBND)
    void FillEstimateSemiAnalytical(const QCluster_t&, Flash_t &) const;

    /// Fills the estimate using the photon library (ICARUS, SBND)
    void FillEstimateLibrary(const QCluster_t&, Flash_t &) const;

  protected:

    void _Configure_(const Config_t &pset);

    double _global_qe;             ///< Global QE for direct light
    double _global_qe_refl;        ///< Global QE for reflected light
    double _sigma_qe;              ///< Sigma for Gaussian centered on Global QE
    std::vector<double> _qe_v;     ///< PMT-wise relative QE
    bool _use_semi_analytical;     ///< If the semi-analytical approach should be used
    #if USING_LARSOFT == 1
    larg4::OpFastScintillation* _opfast_scintillation; ///< For SBND semi-analytical
    #endif
  };

  /**
     \class flashmatch::PhotonLibHypothesisFactory
  */
  class PhotonLibHypothesisFactory : public FlashHypothesisFactoryBase {
  public:
    /// ctor
    PhotonLibHypothesisFactory() { FlashHypothesisFactory::get().add_factory("PhotonLibHypothesis",this); }
    /// dtor
    ~PhotonLibHypothesisFactory() {}
    /// creation method
    BaseFlashHypothesis* create(const std::string instance_name) { return new PhotonLibHypothesis(instance_name); }
  };
}
#endif

/** @} */ // end of doxygen group
