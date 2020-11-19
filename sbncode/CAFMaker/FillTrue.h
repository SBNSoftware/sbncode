
#ifndef CAF_FILLTRUE_H
#define CAF_FILLTRUE_H

#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "TRandom.h"
#include "TDatabasePDG.h"

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"

#include "nusimdata/SimulationBase/GTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCTrack.h"

#include "sbncode/StandardRecord/SRFakeReco.h"
#include "sbncode/StandardRecord/SRTrueParticle.h"
#include "sbncode/StandardRecord/SRTruthMatch.h"
#include "sbncode/StandardRecord/StandardRecord.h"

namespace caf
{


  void FillSliceTruth(const std::vector<art::Ptr<recob::Hit>> &hits,
                      const std::vector<art::Ptr<simb::MCTruth>> &neutrinos,
                      const std::vector<caf::SRTrueInteraction> &srneutrinos,
                      const cheat::ParticleInventoryService &inventory_service,
                      const detinfo::DetectorClocksData &clockData,
                      caf::SRSlice &srslice, caf::SRTruthBranch &srmc,
                      bool allowEmpty = false);

  void FillTrueG4Particle(const simb::MCParticle &mcparticle, 
			  const std::vector<geo::BoxBoundedGeo> &active_volumes,
			  const std::vector<std::vector<geo::BoxBoundedGeo>> &tpc_volumes,
                          const std::map<int, std::vector<const sim::IDE *>> &id_to_ide_map,
			  const cheat::BackTrackerService &backtracker,
			  const cheat::ParticleInventoryService &inventory_service,
			  const std::vector<art::Ptr<simb::MCTruth>> &neutrinos,
                          caf::SRTrueParticle &srparticle); 

  // TODO: implement
  void FillTrueNeutrino(const art::Ptr<simb::MCTruth> mctruth, 
			const simb::MCFlux &mcflux, 
			const std::vector<caf::SRTrueParticle> &srparticles,
			caf::SRTrueInteraction &srneutrino, size_t i);

  void FillTrackTruth(const std::vector<art::Ptr<recob::Hit>> &hits,
                      const std::vector<caf::SRTrueParticle> &particles,
                      const detinfo::DetectorClocksData &clockData,
		      caf::SRTrack& srtrack,
		      bool allowEmpty = false);

  void FillShowerTruth(const std::vector<art::Ptr<recob::Hit>> &hits,
                      const std::vector<caf::SRTrueParticle> &particles,
                      const detinfo::DetectorClocksData &clockData,
		      caf::SRShower& srshower,
		      bool allowEmpty = false);

  void FillFakeReco(const std::vector<art::Ptr<simb::MCTruth>> &mctruths, 
                    const std::vector<art::Ptr<sim::MCTrack>> &mctracks, 
                    const std::vector<geo::BoxBoundedGeo> &volumes,
                    TRandom &rand,
                    std::vector<caf::SRFakeReco> &srfakereco);

  std::map<int, std::vector<const sim::IDE*>> PrepSimChannels(const std::vector<art::Ptr<sim::SimChannel>> &simchannels);

}

#endif
