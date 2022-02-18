#ifndef TRACK_STOPPING_CHI2_ALG_H_SEEN
#define TRACK_STOPPING_CHI2_ALG_H_SEEN

///////////////////////////////////////////////////////////////////////////////
//
// TrackStoppingChi2Alg.h
//
// Algorithm that performs an exponential and a 0-order polynomial fit
// to stopping tracks in order to identify Bragg peaks
//
///////////////////////////////////////////////////////////////////////////////

#include "fhiclcpp/ParameterSet.h"
#include "sbnobj/Common/Reco/StoppingChi2Fit.h"
#include "sbncode/GeometryTools/TPCGeoAlg.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"

#include <vector>

namespace sbn {

  class TrackStoppingChi2Alg {

  public:

    TrackStoppingChi2Alg(fhicl::ParameterSet const& p);

    StoppingChi2Fit RunFit(const std::vector<float> &dEdxVec, const std::vector<float> &resRangeVec) const;
    
    // Prepare dE/dx and residual range vectors for fitting assuming pandora's track direction
    StoppingChi2Fit RunFit(const anab::Calorimetry& calo) const;

    // Prepare dE/dx and residual range vectors for fitting assuming incoming cosmic hypothesis
    StoppingChi2Fit RunFitForCosmicID(const anab::Calorimetry& calo) const;
    
  private:
  
    const float fFitRange, fMaxdEdx;
    const unsigned int fMinHits;

    sbn::TPCGeoAlg fTpcGeo;

  };
}

#endif
