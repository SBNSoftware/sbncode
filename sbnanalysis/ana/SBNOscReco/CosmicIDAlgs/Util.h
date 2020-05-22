#ifndef COSMICIDUTILS_H_SEEN
#define COSMICIDUTILS_H_SEEN


///////////////////////////////////////////////
// CosmicIdUtils.h
//
// Reco utilities for doing cosmic removal in ana modules
// T Brooks (tbrooks@fnal.gov), November 2018
///////////////////////////////////////////////

// LArSoft
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "canvas/Persistency/Common/Ptr.h"

// c++
#include <vector>
#include <utility>
#include <set>

namespace ana{
namespace CosmicIdUtils{

  // Determine if there is a PDS flash in time with the neutrino beam
  bool BeamFlash(std::vector<double> flashes, double beamTimeMin, double beamTimeMax);

  geo::TPCID DetectedInTPC(const std::vector<art::Ptr<recob::Hit>> &hits);
  std::set<geo::TPCID> DetectedInTPCS(const std::vector<art::Ptr<recob::Hit>> &hits);
  
}
}

#endif
