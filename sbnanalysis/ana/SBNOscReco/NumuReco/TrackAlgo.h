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
  void ApplyTrueParticleID(const RecoInteraction &interaction, std::map<size_t, RecoTrack> &tracks, const std::map<size_t, TrueParticle> &particles);
  void ApplyParticleID(const RecoInteraction &interaction, std::map<size_t, RecoTrack> &tracks);
  void ApplyScoredParticleID(const RecoInteraction &interaction, std::map<size_t, RecoTrack> &tracks, float cut);

  float TrackPurity(const numu::RecoTrack &track, const std::map<size_t, numu::TrueParticle> &particles); 
}

#endif
