include_directories(${GENIE_INCLUDE_DIRS})
art_make(
  LIBRARY_NAME sbncode_SBNEventWeight_Calculators_BNBFlux
  LIB_LIBRARIES
    ${GENIE_LIB_LIST}
    sbncode::SBNEventWeight_Base
    nugen::EventGeneratorBase_GENIE
    nugen::NuReweight_art
    nugen::NuReweight
    nurandom::RandomUtils_NuRandomService_service
    art::Framework_Principal
    art::Persistency_Provenance
    art::Framework_Services_Registry
    fhiclcpp::fhiclcpp
    cetlib::cetlib
    nusimdata::SimulationBase
    ROOT::Hist
  MODULE_LIBRARIES
    ${GENIE_LIB_LIST}
)

install_headers()
install_fhicl()
install_source()

