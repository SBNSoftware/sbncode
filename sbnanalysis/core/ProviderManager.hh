#ifndef __sbnanalysis_core_ProviderManager__
#define __sbnanalysis_core_ProviderManager__

#include <string>
#include "core/Experiment.hh"

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
  ProviderManager(Experiment det, std::string fcl="");

  const geo::GeometryCore* GetGeometryProvider() {
    return fGeometryProvider.get();
  }

  const detinfo::LArPropertiesStandard* GetLArPropertiesProvider() {
    return fLArPropertiesProvider.get();
  }

  const detinfo::DetectorClocksStandard* GetDetectorClocksProvider() {
    return fDetectorClocksProvider.get();
  }

  const detinfo::DetectorPropertiesStandard* GetDetectorPropertiesProvider() {
    return fDetectorPropertiesProvider.get();
  }

  const cheat::BackTracker* GetBackTrackerProvider() {
    return fBackTrackerProvider.get();
  }

  const cheat::ParticleInventory* GetParticleInventoryProvider() {
    return fParticleInventoryProvider.get();
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
};

}  // namespace core

#endif  // __sbnanalysis_core_ProviderManager__

