#include <cstdlib>
#include <iostream>
#include <string>
#include "ubcore/Geometry/UBooNEGeometryHelper.h"
#include "larcorealg/Geometry/StandaloneGeometrySetup.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "icaruscode/Geometry/ChannelMapIcarusAlg.h"
#include "sbndcode/Geometry/ChannelMapSBNDAlg.h"
#include "ubcore/Geometry/ChannelMapUBooNEAlg.h"
#include "larcorealg/Geometry/StandaloneBasicSetup.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandardTestHelpers.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardataalg/DetectorInfo/DetectorClocksStandardTestHelpers.h"
#include "lardataalg/DetectorInfo/DetectorClocksStandard.h"
#include "lardataalg/DetectorInfo/LArPropertiesStandardTestHelpers.h"
#include "lardataalg/DetectorInfo/LArPropertiesStandard.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib/filepath_maker.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "core/ServiceManager.hh"

namespace core {

ServiceManager::ServiceManager(Detector det, std::string fcl) {
  // Configuration look up policy: allow absolute paths, and search in the
  // $FHICL_FILE_PATH for non-absolute paths.
  const char* fhicl_file_path = getenv("FHICL_FILE_PATH");
  assert(fhicl_file_path);
  cet::filepath_lookup_nonabsolute policy(fhicl_file_path);

  fhicl::ParameterSet cfg;

  // Geometry: Load the correct geo::GeometryCore subclass for each detector
  std::unique_ptr<geo::GeometryCore> pgeom = NULL;

  switch (det) {
    case kSBND: {
      fcl = (fcl == "" ? "gallery_services_sbnd.fcl" : fcl);
      cfg = lar::standalone::ParseConfiguration(fcl, policy);
      pgeom = \
        lar::standalone::SetupGeometry<geo::ChannelMapSBNDAlg>(
          cfg.get<fhicl::ParameterSet>("services.Geometry"));
      }
      break;
    case kUBOONE: {
      fcl = (fcl == "" ? "gallery_services_uboone.fcl" : fcl);
      cfg = lar::standalone::ParseConfiguration(fcl, policy);
      fhicl::ParameterSet pset = \
        cfg.get<fhicl::ParameterSet>("services.Geometry");
      // For MicroBooNE, we also need to initialize the channel map with the
      // optical detector channel mapping
      std::shared_ptr<geo::ChannelMapUBooNEAlg> channelMap = \
        std::make_shared<geo::ChannelMapUBooNEAlg>(
          cfg.get<fhicl::ParameterSet>("services.ExptGeoHelperInterface"), pset);
      pgeom = \
        lar::standalone::SetupGeometryWithChannelMapping(pset, channelMap);
      }
      break;
    case kICARUS: {
      fcl = (fcl == "" ? "gallery_services_icarus.fcl" : fcl);
      cfg = lar::standalone::ParseConfiguration(fcl, policy);
      pgeom = \
        lar::standalone::SetupGeometry<geo::ChannelMapIcarusAlg>(
          cfg.get<fhicl::ParameterSet>("services.Geometry"));
      }
      break;
    default:
      std::cerr << "ServiceManager: Unknown detector ID" << std::endl;
      assert(false);
      break;
  }

  fGeometryService = pgeom.get();

  // LArProperties
  auto larp = \
    testing::setupProvider<detinfo::LArPropertiesStandard>(
      cfg.get<fhicl::ParameterSet>("services.LArPropertiesService"));
  fLArPropertiesService = larp.get();

  // DetectorClocks
  auto clks = \
    testing::setupProvider<detinfo::DetectorClocksStandard>(
      cfg.get<fhicl::ParameterSet>("services.DetectorClocksService"));
  fDetectorClocksService = clks.get();

  // DetectorProperties
  auto prop = \
    testing::setupProvider<detinfo::DetectorPropertiesStandard>(
      cfg.get<fhicl::ParameterSet>("services.DetectorPropertiesService"),
      detinfo::DetectorPropertiesStandard::providers_type {
        fGeometryService,
        static_cast<detinfo::LArProperties const*>(fLArPropertiesService),
        static_cast<detinfo::DetectorClocks const*>(fDetectorClocksService)
      });
  fDetectorPropertiesService = prop.get();

  config = new fhicl::ParameterSet(cfg);

  std::cout << "ServiceManager: Loaded configuration for: "
            << fGeometryService->DetectorName() << std::endl;
}

}  // namespace core

