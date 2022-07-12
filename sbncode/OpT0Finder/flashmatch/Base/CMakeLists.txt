set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore  -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl)
link_libraries( ${LIB_NAME} -L$ENV{BOOST_LIB} ${ROOTLIB} )
art_make(
    LIB_LIBRARIES
        sbncode_OpT0Finder_flashmatch_GeoAlgo
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
	larsim_PhotonPropagation_PhotonVisibilityService_service
        larsim_Simulation
        lardataalg_DetectorInfo
        art_Framework_Core
        art_Framework_Services_Registry
        art_Framework_Principal
        art_Persistency_Provenance
        art_Utilities
        canvas
        art_Utilities canvas
        ${ART_FRAMEWORK_CORE}
        ${ART_FRAMEWORK_IO_SOURCES}
        ${ART_FRAMEWORK_PRINCIPAL}
        ${ART_FRAMEWORK_SERVICES_REGISTRY}
        ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
        ${ART_ROOT_IO_TFILESERVICE_SERVICE}
        ${ART_PERSISTENCY_COMMON}
        ${ART_PERSISTENCY_PROVENANCE}
        ${ART_UTILITIES}
        ${MF_MESSAGELOGGER}
        ${MF_UTILITIES}
        ${CETLIB}
        ${CLHEP}
        ${BOOST_LIB}
        ${Boost_SYSTEM_LIBRARY}
        ${ROOT_BASIC_LIB_LIST}
        ${FHICLCPP}
        cetlib cetlib_except
)

install_headers()
install_source()
