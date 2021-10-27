#ifndef CAFRECOSELECTIONUTILS_H
#define CAFRECOSELECTIONUTILS_H


///////////////////////////////////////////////
// RecoUtils.h
//
// A few reco utilities like truth matching 
// D Brailsford (adapted from work by D Brailsford and M Wallbank), October 2017
///////////////////////////////////////////////

// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h" 
#include "art/Framework/Principal/Handle.h" 
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
//#include "lardataobj/AnalysisBase/MVAPIDResult.h"
//#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larcore/Geometry/Geometry.h"

// c++
#include <vector>
#include <map>

// ROOT
#include "TTree.h"

namespace CAFRecoUtils{

  std::vector<std::pair<int, float>> AllTrueParticleIDEnergyMatches(const detinfo::DetectorClocksData &clockData, const std::vector<art::Ptr<recob::Hit> >& hits, bool rollup_unsaved_ids=1);
  float TotalHitEnergy(const detinfo::DetectorClocksData &clockData, const std::vector<art::Ptr<recob::Hit> >& hits);

  float TrackPurity(const detinfo::DetectorClocksData &clockData, int mcparticle_id, const std::vector<art::Ptr<recob::Hit>> &reco_track_hits);
  float TrackCompletion(const detinfo::DetectorClocksData &clockData, int mcparticle_id, const std::vector<art::Ptr<recob::Hit>> &reco_track_hits);

}
#endif
