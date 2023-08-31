add_subdirectory(PhysicsLists)

cet_build_plugin(G4InfoReducer art::module
                       LIBRARIES
                           larsim::Simulation
                           larcore::Geometry_Geometry_service
                           lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
                           art::Persistency_Common
                           art::Persistency_Provenance
                           art::Utilities canvas::canvas
                           BASENAME_ONLY)

set(
  MODULE_LIBRARIES
    art::Utilities
    cetlib_except::cetlib_except
    ROOT::Core
    ROOT::Hist
    ROOT::Physics
    ROOT::Geom
    ROOT::MathCore
    CLHEP::CLHEP
    Geant4::G4digits_hits
    Geant4::G4event
    Geant4::G4geometry
    Geant4::G4global
    Geant4::G4graphics_reps
    Geant4::G4materials
    Geant4::G4intercoms
    Geant4::G4particles
    Geant4::G4physicslists
    Geant4::G4processes
    Geant4::G4track
    Geant4::G4run
    sbncode_LArG4_PhysicsLists
)

cet_build_plugin(PhysListLoader art::module SOURCE PhysListLoader_module.cc LIBRARIES ${MODULE_LIBRARIES})

cet_make_library(SOURCE MergeSimSourcesLiteUtility.cxx
  LIBRARIES
  PUBLIC
  lardataobj::MCBase
  larsim::MergeSimSources
)

cet_build_plugin(MergeSimSourcesSBN art::EDProducer
  LIBRARIES
  sbncode::LArG4
  larsim::MergeSimSources
  lardataobj::MCBase
  larcorealg::headers
  art::Framework_Principal
  art::Framework_Services_Registry
  art::Persistency_Common
  messagefacility::MF_MessageLogger
  canvas::canvas
  fhiclcpp::types
  fhiclcpp::fhiclcpp
  BASENAME_ONLY
)

install_headers()
install_fhicl()
install_source()
