#ifndef __sbnanalysis_core_ServiceManager__
#define __sbnanalysis_core_ServiceManager__

#include <string>

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

namespace core {

/** IDs for each SBN detector. */
typedef enum { kSBND, kUBOONE, kICARUS, kNDETS } Detector;

/**
 * \class ServiceManager
 * \brief Interface to LArSoft services
 *
 * The ServiceManager handles loading and provides access to LArSoft services:
 *
 *   * Geometry
 *   * LArProperties
 *   * DetectorClocks
 *   * DetectorProperties
 *
 * \param det The detector ID (as defined in the Detector enum)
 * \param fcl An optional FHiCL file for service settings
 */
class ServiceManager {
public:
  ServiceManager(Detector det, std::string fcl="");

  const geo::GeometryCore* GetGeometryService() {
    return fGeometryService;
  }

  const detinfo::LArPropertiesStandard* GetLArPropertiesService() {
    return fLArPropertiesService;
  }

  const detinfo::DetectorClocksStandard* GetDetectorClocksService() {
    return fDetectorClocksService;
  }

  const detinfo::DetectorPropertiesStandard* GetDetectorPropertiesService() {
    return fDetectorPropertiesService;
  }

private:
  fhicl::ParameterSet* config;
  geo::GeometryCore* fGeometryService;
  detinfo::LArPropertiesStandard* fLArPropertiesService;
  detinfo::DetectorClocksStandard* fDetectorClocksService;
  detinfo::DetectorPropertiesStandard* fDetectorPropertiesService;
};

}  // namespace core

#endif  // __sbnanalysis_core_ServiceManager__

