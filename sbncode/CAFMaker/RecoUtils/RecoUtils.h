#ifndef CAFRECOSELECTIONUTILS_H
#define CAFRECOSELECTIONUTILS_H


///////////////////////////////////////////////
// RecoUtils.h
//
// A few reco utilities like truth matching 
// D Brailsford (adapted from work by D Brailsford and M Wallbank), October 2017
///////////////////////////////////////////////

// framework
#include "canvas/Persistency/Common/Ptr.h" 

// LArSoft
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "lardataobj/RecoBase/Hit.h"

// C++
#include <vector>

namespace CAFRecoUtils{

  std::vector<std::pair<int, float>> AllTrueParticleIDEnergyMatches(const detinfo::DetectorClocksData &clockData, const std::vector<art::Ptr<recob::Hit> >& hits, bool rollup_unsaved_ids=1);
  float TotalHitEnergy(const detinfo::DetectorClocksData &clockData, const std::vector<art::Ptr<recob::Hit> >& hits);

  float TrackPurity(const detinfo::DetectorClocksData &clockData, int mcparticle_id, const std::vector<art::Ptr<recob::Hit>> &reco_track_hits);
  float TrackCompletion(const detinfo::DetectorClocksData &clockData, int mcparticle_id, const std::vector<art::Ptr<recob::Hit>> &reco_track_hits);

  int GetShowerPrimary(const int g4ID);
}
#endif
