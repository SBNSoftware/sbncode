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

#include <iostream>

#ifndef USING_LARSOFT
#define USING_LARSOFT 1
#endif

#if USING_LARSOFT == 0
#include "flashmatch/Base/BaseFlashHypothesis.h"
#include "flashmatch/Base/FlashHypothesisFactory.h"
#include "flashmatch/Base/FMWKInterface.h"
#include "flashmatch/Base/OpT0FinderException.h"
#else
#include "sbncode/OpT0Finder/flashmatch/Base/BaseFlashHypothesis.h"
#include "sbncode/OpT0Finder/flashmatch/Base/FlashHypothesisFactory.h"
#include "sbncode/OpT0Finder/flashmatch/Base/FMWKInterface.h"
#include "sbncode/OpT0Finder/flashmatch/Base/OpT0FinderException.h"
#endif

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

    void BuildHypothesis(const QCluster_t& trk, Flash_t &flash) const;

    int InspectTouchingEdges(const QCluster_t&) const;

    QCluster_t TrackExtension(const QCluster_t&, const int touch) const;

    QCluster_t ComputeExtension(const geoalgo::Vector& A, const geoalgo::Vector& B) const;

    void TrackExtension(const QCluster_t&, Flash_t&) const;

  protected:

    void _Configure_(const Config_t &pset);

    double _global_qe;             ///< Global QE
    double _global_qe_refl;        ///< Global QE for reflected light
    double _sigma_qe;              ///< Sigma for Gaussian centered on Global QE
    std::vector<double> _qe_v;     ///< PMT-wise relative QE
    double _reco_pe_calib;         ///< A global calibration factor for reconstructed PE 
    double _segment_size;
    bool _extend_tracks;
    double _threshold_proximity;
    double _threshold_track_len;
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
