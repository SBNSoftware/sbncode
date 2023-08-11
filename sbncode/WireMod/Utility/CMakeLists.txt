cet_make_library(
  LIBRARY_NAME
    sbncode_WireMod_Utility

  SOURCE
    WireModUtility.cc

  LIBRARIES
    ${ART_FRAMEWORK_CORE}
    ${MF_MESSAGELOGGER}
    ${Boost_SYSTEM_LIBRARY}
    ${ROOT_BASIC_LIB_LIST}
    ${ROOT_HIST}
    art_root_io::TFileService_service
    lardataobj::AnalysisBase
    lardataobj::RawData
    lardataobj::RecoBase
    lardata::RecoObjects
    larreco::RecoAlg
    ${LARRECO_LIB}
    ${LARDATA_LIB}
    ${ART_FRAMEWORK_CORE}
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
    ${MF_MESSAGELOGGER}
    ${FHICLCPP}
    ${CLHEP}
    ${ROOT_GEOM}
    ${ROOT_XMLIO}
    ${ROOT_GDML}
    ${ROOT_BASIC_LIB_LIST}
    ${ROOT_HIST}
    BASENAME_ONLY
)

install_headers()
install_source()

