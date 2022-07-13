art_make(
  LIBRARY_NAME sbncode_SBNEventWeight_Base
  LIB_LIBRARIES
    sbnobj_Common_SBNEventWeight
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
    ${ART_UTILITIES}
    ${ROOT_BASIC_LIB_LIST}
    ${CLHEP}
    canvas
    cetlib_except
)

install_headers()
install_fhicl()
install_source()

