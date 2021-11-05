
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
#include "sbnobj/Common/SBNEventWeight/EventWeightMap.h"
#include "sbnobj/Common/SBNEventWeight/EventWeightParameterSet.h"
#include "lardataobj/MCBase/MCTrack.h"

#include "sbnobj/Common/EventGen/MeVPrtl/MeVPrtlTruth.h"

#include "sbnanaobj/StandardRecord/SRFakeReco.h"
#include "sbnanaobj/StandardRecord/SRGlobal.h"
#include "sbnanaobj/StandardRecord/SRTrueParticle.h"
#include "sbnanaobj/StandardRecord/SRTruthMatch.h"
#include "sbnanaobj/StandardRecord/StandardRecord.h"
#include "sbnanaobj/StandardRecord/SRMeVPrtl.h"

namespace caf
{
  // Helpers
  caf::Wall_t GetWallCross( const geo::BoxBoundedGeo &volume,
        const TVector3 p0,
        const TVector3 p1);

  caf::g4_process_ GetG4ProcessID(const std::string &name);
  
  void FillSRGlobal(const sbn::evwgh::EventWeightParameterSet& pset,
                    caf::SRGlobal& srglobal,
                    std::map<std::string, unsigned int>& weightPSetIndex);

  void FillSliceTruth(const std::vector<art::Ptr<recob::Hit>> &hits,
                      const std::vector<art::Ptr<simb::MCTruth>> &neutrinos,
                      const caf::SRTruthBranch &srmc,
                      const cheat::ParticleInventoryService &inventory_service,
                      const detinfo::DetectorClocksData &clockData,
                      caf::SRSlice &srslice, 
                      bool allowEmpty = false);

  void FillSliceFakeReco(const std::vector<art::Ptr<recob::Hit>> &hits,
                         const std::vector<art::Ptr<simb::MCTruth>> &neutrinos,
                         const caf::SRTruthBranch &srmc,
                         const cheat::ParticleInventoryService &inventory_service,
                         const detinfo::DetectorClocksData &clockData,
                         caf::SRSlice &srslice, 
                         const std::vector<art::Ptr<sim::MCTrack>> &mctracks,
                         const std::vector<geo::BoxBoundedGeo> &volumes, TRandom &rand);

  void FillTrueG4Particle(const simb::MCParticle &particle,
        const std::vector<geo::BoxBoundedGeo> &active_volumes,
        const std::vector<std::vector<geo::BoxBoundedGeo>> &tpc_volumes,
        const std::map<int, std::vector<std::pair<geo::WireID, const sim::IDE *>>> &id_to_ide_map,
        const std::map<int, std::vector<art::Ptr<recob::Hit>>> &id_to_truehit_map,
        const cheat::BackTrackerService &backtracker,
        const cheat::ParticleInventoryService &inventory_service,
        const std::vector<art::Ptr<simb::MCTruth>> &neutrinos,
                          caf::SRTrueParticle &srparticle);

  void FillMeVPrtlTruth(const evgen::ldm::MeVPrtlTruth &truth,
                        const std::vector<geo::BoxBoundedGeo> &active_volumes,
                        caf::SRMeVPrtl &srtruth);

  void FillTrueNeutrino(const art::Ptr<simb::MCTruth> mctruth, 
			const simb::MCFlux &mcflux, 
                        const simb::GTruth& gtruth,
			const std::vector<caf::SRTrueParticle> &srparticles,
                        const std::map<int, std::vector<art::Ptr<recob::Hit>>> &id_to_truehit_map,
			caf::SRTrueInteraction &srneutrino, size_t i,
                        const std::vector<geo::BoxBoundedGeo> &active_volumes);

  void FillEventWeight(const sbn::evwgh::EventWeightMap& wgtmap,
                       caf::SRTrueInteraction& srint,
                       const std::map<std::string, unsigned int>& weightPSetIndex);

  void FillTrackTruth(const std::vector<art::Ptr<recob::Hit>> &hits,
                      const std::vector<caf::SRTrueParticle> &particles,
                      const detinfo::DetectorClocksData &clockData,
		      caf::SRTrack& srtrack,
		      bool allowEmpty = false);

  void FillStubTruth(const std::vector<art::Ptr<recob::Hit>> &hits,
                     const std::vector<caf::SRTrueParticle> &particles,
                     const detinfo::DetectorClocksData &clockData,
                     caf::SRStub& srstub,
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

  std::map<int, std::vector<std::pair<geo::WireID, const sim::IDE*>>> PrepSimChannels(const std::vector<art::Ptr<sim::SimChannel>> &simchannels, const geo::GeometryCore &geo);
  std::map<int, std::vector<art::Ptr<recob::Hit>>> PrepTrueHits(const std::vector<art::Ptr<recob::Hit>> &allHits, 
    const detinfo::DetectorClocksData &clockData, const cheat::BackTrackerService &backtracker);

}

#endif
