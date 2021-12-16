art_make(
  LIBRARY_NAME sbncode_SBNEventWeight_Base
  LIB_LIBRARIES
    sbnobj_Common_SBNEventWeight
    art::Framework_Services_Registry
    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
    art::Framework_Services_Optional_RandomNumberGenerator_service
    art::Utilities
    ${ROOT_BASIC_LIB_LIST}
    CLHEP::CLHEP
    canvas::canvas
    cetlib_except::cetlib_except
)

install_headers()
install_fhicl()
install_source()

