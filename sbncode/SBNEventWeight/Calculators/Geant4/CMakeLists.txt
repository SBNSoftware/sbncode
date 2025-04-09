include_directories ( $ENV{GEANT4_FQ_DIR}/include )
include_directories ( $ENV{GEANT4REWEIGHT_INC} )

art_make(
  LIBRARY_NAME sbncode_SBNEventWeight_Calculators_Geant4
  LIBRARIES
    sbncode_SBNEventWeight_Base
    sbncode_SBNEventWeight_Calculators_MultiVariateRNG
    nusimdata::SimulationBase
    nurandom::RandomUtils_NuRandomService_service
    art::Framework_Principal
    art::Framework_Services_Registry
    art_root_io::TFileService_service
    larcore::Geometry_Geometry_service
    cetlib_except::cetlib_except
    geant4reweight::ReweightBaseLib
    geant4reweight::PropBaseLib
    ROOT::Tree
)

install_headers()
install_fhicl()
install_source()

