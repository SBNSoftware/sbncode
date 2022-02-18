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
#include "lardataobj/AnalysisBase/Calorimetry.h"

#include <vector>

namespace sbn {

  class TrackStoppingChi2Alg {

  public:

    TrackStoppingChi2Alg(fhicl::ParameterSet const& p);

    StoppingChi2Fit RunFit(const std::vector<float> &dEdxVec, const std::vector<float> &resRangeVec) const;
    
    StoppingChi2Fit RunFit(const anab::Calorimetry& calo) const;
    
  private:
  
    const float fFitRange, fMaxdEdx;
    const unsigned int fMinHits;

  };
}

#endif
