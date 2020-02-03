#ifndef _sbncode_TrackAlgo_hh
#define _sbncode_TrackAlgo_hh

#include <map>

#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "../Data/RecoTrack.h"
#include "../Data/RecoEvent.h"

namespace numu {
  float MeanTruncateddQdx(const anab::Calorimetry &calo);
  float TrackMomentum(const numu::RecoTrack &track);
  float RangeMomentum(const numu::RecoTrack &track);
  float MCSMomentum(const numu::RecoTrack &track);
  void ApplyParticleID(const std::map<size_t, numu::RecoTrack> &tracks, const numu::RecoInteraction &reco); 
  void ApplyTrueParticleID(const RecoSlice &slice, std::map<size_t, RecoTrack> &tracks, const std::map<size_t, TrueParticle> &particles);
  void ApplyParticleID(const RecoSlice &slice, std::map<size_t, RecoTrack> &tracks);
}

#endif
