include_directories ( $ENV{GEANT4_FQ_DIR}/include )
include_directories ( $ENV{GEANT4REWEIGHT_INC} )

art_make_library(
  LIBRARY_NAME sbncode_SBNEventWeight_Calculators_Geant4
  LIBRARIES
    sbncode_SBNEventWeight_Base
    nusimdata::SimulationBase
    nurandom::RandomUtils_NuRandomService_service
    art::Framework_Principal
    art::Framework_Services_Registry
    art_root_io::TFileService_service
    larcore::Geometry_Geometry_service
    cetlib_except::cetlib_except
    ReweightBaseLib
    PropBaseLib
    ROOT::Tree
)

install_headers()
install_fhicl()
install_source()

