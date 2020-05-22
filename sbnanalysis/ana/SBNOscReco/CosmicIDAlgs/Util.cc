#include "Util.h"

namespace ana {

// =============================== UTILITY FUNCTIONS ==============================

  // Determine if there is a PDS flash in time with the neutrino beam
  bool CosmicIdUtils::BeamFlash(std::vector<double> flashes, double beamTimeMin, double beamTimeMax){
    //
    bool beamFlash = false;
    std::sort(flashes.begin(), flashes.end());
    // Loop over flashes in tpc 0
    for(size_t i = 0; i < flashes.size(); i++){
      double time = flashes[i];
      if(time > beamTimeMin && time < beamTimeMax) beamFlash = true;
    }

    return beamFlash;
  }

  std::set<geo::TPCID> CosmicIdUtils::DetectedInTPCS(const std::vector<art::Ptr<recob::Hit>> &hits) {
   std::set<geo::TPCID> ret; 
    for (const art::Ptr<recob::Hit> hit: hits) {
      ret.insert(hit->WireID().asTPCID());
    }
    return ret;
  }


  geo::TPCID CosmicIdUtils::DetectedInTPC(const std::vector<art::Ptr<recob::Hit>> &hits) {
    geo::TPCID ret; // constructor defaults to invalid
    for (const art::Ptr<recob::Hit> hit: hits) {
      if (!ret) ret = hit->WireID().asTPCID();
      else if (ret != hit->WireID().asTPCID()) {
        geo::TPCID invalid;
        return invalid;
      }
    }
    return ret;
  }

}
