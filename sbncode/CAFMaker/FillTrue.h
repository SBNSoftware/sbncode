
#ifndef CAF_FILLTRUE_H
#define CAF_FILLTRUE_H

#include "art/Framework/Services/Registry/ServiceHandle.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "nusimdata/SimulationBase/MCParticle.h"

#include "sbncode/StandardRecord/SRTrueParticle.h"
#include "sbncode/StandardRecord/StandardRecord.h"

namespace caf
{
  struct ParticleData {
    const std::vector<geo::BoxBoundedGeo> *AV;
    const std::vector<std::vector<geo::BoxBoundedGeo>> *TPCVolumes;
    const cheat::BackTrackerService *backtracker;
    const cheat::ParticleInventoryService *inventory_service;
    std::vector<art::Ptr<simb::MCTruth>> neutrinos;
  };

  void FillTrueG4Particle(const simb::MCParticle &mcparticle, 
                          const ParticleData &pdata,
                          caf::SRTrueParticle &srparticle); 

}

#endif
