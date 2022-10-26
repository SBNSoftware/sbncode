
art_make( 
  LIBRARY_NAME sbncode_SBNEventWeight_Calculators_CrossSection
  LIB_LIBRARIES
    sbncode::SBNEventWeight_Base
    nugen::EventGeneratorBase_GENIE
    nugen::NuReweight_art
    nugen::NuReweight
    nusimdata::SimulationBase
    nurandom::RandomUtils_NuRandomService_service
    art::Framework_Principal
    art::Framework_Services_Registry
    ${GENIE_LIB_LIST}
    log4cpp::log4cpp
    GSL::gsl  # FWNUM relies on GSL
    cetlib_except::cetlib_except
)

install_headers()
install_fhicl()
install_source()

