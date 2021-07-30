#ifndef CAF_FILLTRUE_H
#define CAF_FILLTRUE_H

#include "canvas/Persistency/Common/Ptr.h"

#include <map>
#include <string>
#include <vector>

class TRandom;

namespace cheat
{
  class BackTrackerService;
  class ParticleInventoryService;
}

namespace detinfo
{
  class DetectorClocksData;
}

namespace geo
{
  class BoxBoundedGeo;
  class GeometryCore;
  class WireID;
}

namespace simb
{
  class GTruth;
  class MCFlux;
  class MCNeutrino;
  class MCParticle;
  class MCTruth;
}

namespace sim
{
  class IDE;
  class MCTrack;
  class SimChannel;
}

namespace sbn::evwgh
{
  typedef std::map<std::string, std::vector<float> > EventWeightMap;
  class EventWeightParameterSet;
}

namespace recob
{
  class Hit;
}

namespace caf
{
  class SRFakeReco;
  class SRGlobal;
  class SRShower;
  class SRSlice;
  class SRTrack;
  class SRTrueParticle;
  class SRTruthBranch;
  class SRTruthMatch;
  class SRTrueInteraction;
}

namespace caf
{
  void FillSRGlobal(const sbn::evwgh::EventWeightParameterSet& pset,
                    caf::SRGlobal& srglobal,
                    std::map<std::string, unsigned int>& weightPSetIndex);

  void FillSliceTruth(const std::vector<art::Ptr<recob::Hit>> &hits,
                      const std::vector<art::Ptr<simb::MCTruth>> &neutrinos,
                      const std::vector<caf::SRTrueInteraction> &srneutrinos,
                      const cheat::ParticleInventoryService &inventory_service,
                      const detinfo::DetectorClocksData &clockData,
                      caf::SRSlice &srslice, caf::SRTruthBranch &srmc,
                      bool allowEmpty = false);

  void FillSliceFakeReco(const std::vector<art::Ptr<recob::Hit>> &hits,
                         const std::vector<art::Ptr<simb::MCTruth>> &neutrinos,
                         const std::vector<caf::SRTrueInteraction> &srneutrinos,
                         const cheat::ParticleInventoryService &inventory_service,
                         const detinfo::DetectorClocksData &clockData,
                         caf::SRSlice &srslice, caf::SRTruthBranch &srmc,
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

  void FillTrueNeutrino(const art::Ptr<simb::MCTruth> mctruth, 
			const simb::MCFlux &mcflux, 
                        const simb::GTruth& gtruth,
			const std::vector<caf::SRTrueParticle> &srparticles,
                        const std::map<int, std::vector<art::Ptr<recob::Hit>>> &id_to_truehit_map,
			caf::SRTrueInteraction &srneutrino, size_t i);

  void FillEventWeight(const sbn::evwgh::EventWeightMap& wgtmap,
                       caf::SRTrueInteraction& srint,
                       const std::map<std::string, unsigned int>& weightPSetIndex);

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

  std::map<int, std::vector<std::pair<geo::WireID, const sim::IDE*>>> PrepSimChannels(const std::vector<art::Ptr<sim::SimChannel>> &simchannels, const geo::GeometryCore &geo);
  std::map<int, std::vector<art::Ptr<recob::Hit>>> PrepTrueHits(const std::vector<art::Ptr<recob::Hit>> &allHits, 
    const detinfo::DetectorClocksData &clockData, const cheat::BackTrackerService &backtracker);

}

#endif
