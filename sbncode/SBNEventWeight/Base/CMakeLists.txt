art_make_library(
  LIBRARY_NAME sbncode_SBNEventWeight_Base
  LIBRARIES
    sbnobj::Common_SBNEventWeight
    nurandom::RandomUtils_NuRandomService_service
    art::Framework_Services_Registry
    lardataobj::Simulation
    CLHEP::CLHEP
    canvas::canvas
    cetlib_except::cetlib_except
)
install_headers()
install_fhicl()
install_source()

