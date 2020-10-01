#ifndef _sbncode_TrackAlgo_hh
#define _sbncode_TrackAlgo_hh

#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "../Data/RecoTrack.h"

namespace numu {
  float MeanTruncateddQdx(const anab::Calorimetry &calo);
  float TrackMomentum(const numu::RecoTrack &track);
  float RangeMomentum(const numu::RecoTrack &track);
  float MCSMomentum(const numu::RecoTrack &track);
}

#endif
