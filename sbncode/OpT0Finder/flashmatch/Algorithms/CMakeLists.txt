set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMinuit -pthread -lm -ldl)

link_libraries( ${LIB_NAME} -L$ENV{BOOST_LIB} -lboost_system ${ROOTLIB} )

art_make(
    LIB_LIBRARIES
        larsim_PhotonPropagation
        larsim_PhotonPropagation_PhotonVisibilityService_service
        larsim_LegacyLArG4
        larcorealg_GeoAlgo
        # sbncode_OpT0Finder_flashmatch_GeoAlgo
        sbncode_OpT0Finder_flashmatch_Base
        ${ART_FRAMEWORK_CORE}
        ${ART_FRAMEWORK_PRINCIPAL}
        ${ART_FRAMEWORK_BASIC}
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
        ${ART_PERSISTENCY_COMMON}
        ${ART_PERSISTENCY_PROVENANCE}
        ${ART_UTILITIES}
        ${ROOT_BASIC_LIB_LIST}
        ${FHICLCPP}
        ${CETLIB}
)

install_headers()
install_source()
