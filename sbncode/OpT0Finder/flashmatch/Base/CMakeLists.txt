set( ROOTLIB -L$ENV{ROOTSYS}/lib -lCore  -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl)
link_libraries( ${LIB_NAME} -L$ENV{BOOST_LIB} ${ROOTLIB} )
art_make(
    LIB_LIBRARIES
        larcorealg_GeoAlgo
        # sbncode_OpT0Finder_flashmatch_GeoAlgo
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
        lardataalg_DetectorInfo
        art_Framework_Core
        art::Framework_Services_Registry
        art::Framework_Principal
        art::Persistency_Provenance
        art::Utilities
        canvas::canvas
        art::Utilities canvas::canvas
        art::Framework_Core
        art::Framework_IO_Sources
        art::Framework_Principal
        art::Framework_Services_Registry
        art_root_io::tfile_support ROOT::Core
        art_root_io::TFileService_service
        art::Persistency_Common
        art::Persistency_Provenance
        art::Utilities
        messagefacility::MF_MessageLogger
        
        cetlib::cetlib
        CLHEP::CLHEP
        ${BOOST_LIB}
        
        ${ROOT_BASIC_LIB_LIST}
        fhiclcpp::fhiclcpp
        cetlib::cetlib cetlib_except
)

install_headers()
install_source()
