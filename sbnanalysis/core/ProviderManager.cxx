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
#include "larsim/MCCheater/BackTracker.h"
#include "larsim/MCCheater/ParticleInventory.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib/filepath_maker.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "core/ProviderManager.hh"
#include "core/Experiment.hh"

namespace core {

ProviderManager::ProviderManager(Experiment det, std::string fcl) {
  // Configuration look up policy: allow absolute paths, and search in the
  // $FHICL_FILE_PATH for non-absolute paths.
  const char* fhicl_file_path = getenv("FHICL_FILE_PATH");
  assert(fhicl_file_path);
  cet::filepath_lookup_nonabsolute policy(fhicl_file_path);

  fhicl::ParameterSet cfg;

  // Geometry: Load the correct geo::GeometryCore subclass for each detector
  switch (det) {
    case kExpSBND: {
      fcl = (fcl == "" ? "gallery_services_sbnd.fcl" : fcl);
      cfg = lar::standalone::ParseConfiguration(fcl, policy);
      fGeometryProvider = \
        lar::standalone::SetupGeometry<geo::ChannelMapSBNDAlg>(
          cfg.get<fhicl::ParameterSet>("services.Geometry"));
      }
      break;
    case kExpMicroBooNE: {
      fcl = (fcl == "" ? "gallery_services_uboone.fcl" : fcl);
      cfg = lar::standalone::ParseConfiguration(fcl, policy);
      fhicl::ParameterSet pset = \
        cfg.get<fhicl::ParameterSet>("services.Geometry");
      // For MicroBooNE, we also need to initialize the channel map with the
      // optical detector channel mapping
      std::shared_ptr<geo::ChannelMapUBooNEAlg> channelMap = \
        std::make_shared<geo::ChannelMapUBooNEAlg>(
          cfg.get<fhicl::ParameterSet>("services.ExptGeoHelperInterface"), pset);
      fGeometryProvider = \
        lar::standalone::SetupGeometryWithChannelMapping(pset, channelMap);
      }
      break;
    case kExpICARUS: {
      fcl = (fcl == "" ? "gallery_services_icarus.fcl" : fcl);
      cfg = lar::standalone::ParseConfiguration(fcl, policy);
      fGeometryProvider = \
        lar::standalone::SetupGeometry<geo::ChannelMapIcarusAlg>(
          cfg.get<fhicl::ParameterSet>("services.Geometry"));
      }
      break;
    default:
      std::cerr << "ProviderManager: Unknown detector ID" << std::endl;
      assert(false);
      break;
  }

  // LArProperties
 fLArPropertiesProvider = \
    testing::setupProvider<detinfo::LArPropertiesStandard>(
      cfg.get<fhicl::ParameterSet>("services.LArPropertiesService"));

  // DetectorClocks
  fDetectorClocksProvider = \
    testing::setupProvider<detinfo::DetectorClocksStandard>(
      cfg.get<fhicl::ParameterSet>("services.DetectorClocksService"));

  // DetectorProperties
  fDetectorPropertiesProvider = \
    testing::setupProvider<detinfo::DetectorPropertiesStandard>(
      cfg.get<fhicl::ParameterSet>("services.DetectorPropertiesService"),
      detinfo::DetectorPropertiesStandard::providers_type {
        fGeometryProvider.get(),
        static_cast<const detinfo::LArProperties*>(fLArPropertiesProvider.get()),
        static_cast<const detinfo::DetectorClocks*>(fDetectorClocksProvider.get())
      });

  // ParticleInventory
  fParticleInventoryProvider = NULL;
  if (cfg.has_key("services.ParticleInventoryService")) {
    fParticleInventoryProvider = testing::setupProvider<cheat::ParticleInventory>(
       cfg.get<fhicl::ParameterSet>("services.ParticleInventoryService"));
  }
  else {
    std::cerr << "Warning: Particle inventory service is missing from fhicl config (" << fcl << ")." \
      << " Setting ParticleInventoryService to NULL" << std::endl;
  }

  // BackTracker
  fBackTrackerProvider = NULL;
  if (cfg.has_key("services.BackTrackerService")) {
    fBackTrackerProvider = \
      testing::setupProvider<cheat::BackTracker>(
        cfg.get<fhicl::ParameterSet>("services.BackTrackerService"),
        fParticleInventoryProvider.get(), fGeometryProvider.get(),
        fDetectorClocksProvider.get());
  }
  else {
    std::cerr << "Warning: BackTracker service is missing from fhicl config (" << fcl << ")." \
      << " Setting BackTrackerService to NULL" << std::endl;

  }

  config = new fhicl::ParameterSet(cfg);

  std::cout << "ProviderManager: Loaded configuration for: "
            << fGeometryProvider.get()->DetectorName() << std::endl;
}


std::vector<Experiment> ProviderManager::GetValidExperiments() {
  static std::vector<Experiment> ex = {
    kExpSBND, kExpMicroBooNE, kExpICARUS
  };
  return ex;
}

}  // namespace core

