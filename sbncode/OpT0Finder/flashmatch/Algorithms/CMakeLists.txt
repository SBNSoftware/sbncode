set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMinuit -pthread -lm -ldl)

link_libraries( ${LIB_NAME} -L$ENV{BOOST_LIB} -lboost_system ${ROOTLIB} )

art_make(
    LIB_LIBRARIES
        larsim_PhotonPropagation
        larsim_PhotonPropagation_PhotonVisibilityService_service
        larsim_LegacyLArG4
        sbncode_OpT0Finder_flashmatch_GeoAlgo
        sbncode_OpT0Finder_flashmatch_Base
        art::Framework_Core
        art::Framework_Principal
        ${ART_FRAMEWORK_BASIC}
        art::Framework_Services_Registry
        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
        
        art::Persistency_Common
        art::Persistency_Provenance
        art::Utilities
        ${ROOT_BASIC_LIB_LIST}
        fhiclcpp::fhiclcpp
        cetlib::cetlib
)

install_headers()
install_source()
