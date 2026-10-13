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
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "canvas/Persistency/Common/Ptr.h" 
#include "canvas/Persistency/Common/PtrVector.h" 
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "sbnanaobj/StandardRecord/SREnums.h"
//#include "lardataobj/AnalysisBase/MVAPIDResult.h"
//#include "lardataobj/AnalysisBase/ParticleID.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/WireReadoutGeom.h"

// c++
#include <vector>
#include <map>

// ROOT
#include "TTree.h"

namespace sim {
  class IDE;
}

namespace sbn {
  struct ReadoutIDE {
    geo::WireID wire;              ///< Wire on a given plane closest to the drift path of the charge.
    unsigned short tick = 0;       ///< Time tick at which the charge passes closest to the wire.
    const sim::IDE *ide = nullptr; ///< Deposited charge information.
  };
  
  // Helpers
  std::map<int, std::vector<sbn::ReadoutIDE>> PrepSimChannels(const std::vector<art::Ptr<sim::SimChannel>> &simchannels, const geo::WireReadoutGeom &wireReadout);
  std::map<int, std::vector<art::Ptr<recob::Hit>>> PrepTrueHits(const std::vector<art::Ptr<recob::Hit>> &allHits,
    const detinfo::DetectorClocksData &clockData, const cheat::BackTrackerService &backtracker);
  caf::Wall_t GetWallCross( const geo::BoxBoundedGeo &volume, const TVector3 p0, const TVector3 p1);
  caf::g4_process_ GetG4ProcessID(const std::string &name);
}

namespace CAFRecoUtils{

  std::vector<std::pair<int, float>> AllTrueParticleIDEnergyMatches(const detinfo::DetectorClocksData &clockData, const std::vector<art::Ptr<recob::Hit> >& hits, bool rollup_unsaved_ids=1);
  float TotalHitEnergy(const detinfo::DetectorClocksData &clockData, const std::vector<art::Ptr<recob::Hit> >& hits);

  float TrackPurity(const detinfo::DetectorClocksData &clockData, int mcparticle_id, const std::vector<art::Ptr<recob::Hit>> &reco_track_hits);
  float TrackCompletion(const detinfo::DetectorClocksData &clockData, int mcparticle_id, const std::vector<art::Ptr<recob::Hit>> &reco_track_hits);

  int GetShowerPrimary(const int g4ID);
}



#endif
