cet_build_plugin(DeletePi0Photon art::EDFilter
  LIBRARIES PRIVATE
  larcore::Geometry_Geometry_service
  larcoreobj::SummaryData
  lardataobj::Simulation
  nurandom::RandomUtils_NuRandomService_service
  nusimdata::SimulationBase
  dk2nu::Tree
  art::Framework_Principal
  art::Framework_Services_Registry
  art::Persistency_Common
  canvas::canvas
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  cetlib_except::cetlib_except
  CLHEP::Random
  ROOT::Physics
)

install_headers()
install_fhicl()
install_source()
