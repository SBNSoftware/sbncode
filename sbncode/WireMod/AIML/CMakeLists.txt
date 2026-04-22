cet_make_library(
  LIBRARY_NAME
    sbncode_WireMod_AIML

  SOURCE
    WireModInfer.cc

  LIBRARIES
    ${ART_FRAMEWORK_CORE}
    ${MF_MESSAGELOGGER}
    ${Boost_SYSTEM_LIBRARY}
    ${ROOT_BASIC_LIB_LIST}
    ${ROOT_HIST}
    larcorealg::Geometry
    larcore::Geometry_Geometry_service
    lardataobj::AnalysisBase
    lardataobj::RawData
    lardataobj::RecoBase
    lardataalg::DetectorInfo 
    lardata::RecoObjects
    lardata::Utilities
    larreco::RecoAlg
    larsim::MCCheater_ParticleInventoryService_service
    larsim::Utils
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
    Eigen3::Eigen
    BASENAME_ONLY
)

#set_property(SOURCE WireModInfer.cc APPEND PROPERTY COMPILE_DEFINITIONS EIGEN_FFTW_DEFAULT)

install_headers()
install_source()

