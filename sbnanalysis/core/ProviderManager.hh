#ifndef __sbnanalysis_core_ProviderManager__
#define __sbnanalysis_core_ProviderManager__

#include <string>
#include <memory>
#include <vector>
#include "core/Experiment.hh"

#include "gallery/Event.h"

namespace geo {
  class GeometryCore;
}

namespace detinfo {
  class LArPropertiesStandard;
  class DetectorClocksStandard;
  class DetectorPropertiesStandard;
}

namespace fhicl {
  class ParameterSet;
}

namespace cheat {
  class BackTracker;
  class ParticleInventory;
  class PhotonBackTracker;
}

namespace core {

/**
 * \class ProviderManager
 * \brief Interface to LArSoft services
 *
 * The ProviderManager handles loading and provides access to LArSoft services:
 *
 *   * Geometry
 *   * LArProperties
 *   * DetectorClocks
 *   * DetectorProperties
 *
 * \param det The detector ID (as defined in the Detector enum)
 * \param fcl An optional FHiCL file for service settings
 */
class ProviderManager {
public:
  ProviderManager(Experiment det, std::string fcl="", bool setup_event_services=true);

  void SetupServices(gallery::Event &ev);

  const geo::GeometryCore* GetGeometryProvider() const {
    return fGeometryProvider.get();
  }

  const detinfo::LArPropertiesStandard* GetLArPropertiesProvider() const {
    return fLArPropertiesProvider.get();
  }

  const detinfo::DetectorClocksStandard* GetDetectorClocksProvider() const {
    return fDetectorClocksProvider.get();
  }

  const detinfo::DetectorPropertiesStandard* GetDetectorPropertiesProvider() const {
    return fDetectorPropertiesProvider.get();
  }

  cheat::BackTracker* GetBackTrackerProvider() const {
    return fBackTrackerProvider.get();
  }

  cheat::ParticleInventory* GetParticleInventoryProvider() const {
    return fParticleInventoryProvider.get();
  }

  cheat::PhotonBackTracker* GetPhotonBackTrackerProvider() const {
    return fPhotonBackTrackerProvider.get();
  }

  static std::vector<Experiment> GetValidExperiments();

private:
  fhicl::ParameterSet* config;
  std::unique_ptr<geo::GeometryCore> fGeometryProvider;
  std::unique_ptr<detinfo::LArPropertiesStandard> fLArPropertiesProvider;
  std::unique_ptr<detinfo::DetectorClocksStandard> fDetectorClocksProvider;
  std::unique_ptr<detinfo::DetectorPropertiesStandard> fDetectorPropertiesProvider;
  std::unique_ptr<cheat::BackTracker> fBackTrackerProvider;
  std::unique_ptr<cheat::ParticleInventory> fParticleInventoryProvider;
  std::unique_ptr<cheat::PhotonBackTracker> fPhotonBackTrackerProvider;
};

}  // namespace core

#endif  // __sbnanalysis_core_ProviderManager__

