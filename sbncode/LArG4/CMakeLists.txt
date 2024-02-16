cet_build_plugin(G4InfoReducer art::module
                       LIBRARIES
                           larsim::Simulation
                           larcore::Geometry_Geometry_service
                           lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
                           art::Persistency_Common
                           art::Persistency_Provenance
                           art::Utilities canvas::canvas
                           BASENAME_ONLY)

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

cet_build_plugin(SimpleMerge art::EDProducer
  LIBRARIES
  lardataobj::MCBase
  larcorealg::headers
  messagefacility::MF_MessageLogger
)

install_headers()
install_fhicl()
install_source()
