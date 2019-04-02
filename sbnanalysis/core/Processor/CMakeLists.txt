art_make(
  LIBRARY_NAME sbncode_sbnanalysis_core_Base_Processor
  LIB_LIBRARIES
    sbncode_sbnanalysis_Util
    ${ART_FRAMEWORK_CORE}
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
    nusimdata_SimulationBase
    art_Persistency_Common
    art_Persistency_Provenance
    art_Utilities
    canvas
    cetlib cetlib_except
    ${FHICLCPP}
    ${CLHEP}
    ${ROOT_BASIC_LIB_LIST}
    ${Boost_SYSTEM_LIBRARY}
)

install_headers()
install_fhicl()
install_source()

