add_subdirectory(DataTypes)
add_subdirectory(fcl)

art_make(
  MODULE_LIBRARIES 
    nusimdata_SimulationBase
    ${ART_FRAMEWORK_CORE}
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_BASIC}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
    art_Persistency_Common canvas
    art_Persistency_Provenance canvas
    art_Utilities canvas
    ${MF_MESSAGELOGGER}
    ${MF_UTILITIES}
    ${FHICLCPP}
    cetlib cetlib_except
    ${ROOT_BASIC_LIB_LIST}
)

install_headers()
install_fhicl()
install_source()

