art_make(
  LIBRARY_NAME sbncode_sbnanalysis_ana_ExampleAnalysis
  LIB_LIBRARIES
    sbnanalysis_core_DataTypes_dict
    sbncode_sbnanalysis_core_Base_Processor
    sbncode_sbnanalysis_Util
    larsim_EventWeight_Base
    lardata_Utilities
    nutools_NuReweight_art
    nutools_NuReweight
    larsim_Simulation
    nutools_RandomUtils_NuRandomService_service
    nusimdata_SimulationBase
    ${ART_FRAMEWORK_CORE}
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
    canvas
    art_Persistency_Common
    art_Persistency_Provenance
    art_Utilities
    cetlib cetlib_except
    ${FHICLCPP}
    ${CLHEP}
    ${ROOT_BASIC_LIB_LIST}
    ${Boost_SYSTEM_LIBRARY}
)

install_headers()
install_fhicl()
install_source()

