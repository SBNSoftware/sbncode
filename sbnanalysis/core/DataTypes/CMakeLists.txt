include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(ROOTLIB -L$ENV{ROOTSYS}/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl)

link_libraries( ${LIB_NAME} -L$ENV{BOOST_LIB} -lboost_system ${ROOTLIB} )

art_make(
  LIBRARY_NAME sbncode_sbnanalysis_core_DataTypes
  LIB_LIBRARIES
    lardataobj_RawData
    lardataobj_RecoBase
    lardataobj_AnalysisBase
    lardata_Utilities
    ${LARDATA_LIB}
    ${ART_FRAMEWORK_CORE}
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_FRAMEWORK_BASIC}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL}
    ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
    art_Persistency_Common canvas_Persistency_Common
    art_Persistency_Provenance canvas_Persistency_Provenance
    art_Utilities canvas_Utilities
    ${ART_PERSISTENCY_COMMON}
    ${ART_PERSISTENCY_PROVENANCE}
    ${ART_UTILITIES}
    ${MF_MESSAGELOGGER}
    ${MF_UTILITIES}
    ${FHICLCPP}
    ${CETLIB}
    ${ROOT_BASIC_LIB_LIST}
)

install_headers()
install_fhicl()
install_source()

