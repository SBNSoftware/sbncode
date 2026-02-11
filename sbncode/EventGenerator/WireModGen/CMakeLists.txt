cet_build_plugin(WireModMuon art::EDProducer
  SOURCE CORSIKAGenSBN.cc WireModMuon_module.cc
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  larcoreobj::SummaryData
  nurandom::RandomUtils_NuRandomService_service
  nusimdata::SimulationBase
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ifdh::ifdh
  messagefacility::MF_MessageLogger
  ROOT::Core
  ROOT::EG
  CLHEP::Random
  SQLite::SQLite3
  lardata::DetectorInfoServices_DetectorPropertiesServiceStandard_service
  lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
  lardata::DetectorClocksService
)

cet_build_plugin(FilterPrimarySimChannel art::EDProducer
  LIBRARIES PRIVATE
  larsim::Simulation
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  messagefacility::MF_MessageLogger
)

cet_build_plugin(FilterPrimarySimEnergyDepo art::EDProducer
  LIBRARIES PRIVATE
  art::Framework_Services_Registry
  larsim::Simulation
  art::Framework_Principal
  art::Framework_Services_Registry
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  messagefacility::MF_MessageLogger
  larcore::Geometry_Geometry_service
  larcorealg::Geometry
  lardata::DetectorInfoServices_DetectorPropertiesServiceStandard_service
  lardata::DetectorInfoServices_DetectorClocksServiceStandard_service
  lardata::DetectorClocksService
)

install_headers()
install_fhicl()
install_source()
