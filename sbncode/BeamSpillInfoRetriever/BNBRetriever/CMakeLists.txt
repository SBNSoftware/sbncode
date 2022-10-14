find_package(ifbeam)
find_package(ifbeam_art)

include_directories( $ENV{IFBEAM_DIR}/include )
include_directories( $ENV{LIBWDA_DIR}/include )
link_directories( $ENV{LIBWDA_DIR}/lib )
include_directories( $ENV{IFDHC_DIR}/inc )
link_directories( $ENV{IFDHC_DIR}/lib )
include_directories( $ENV{IFDH_ART_DIR}/include )
include_directories( $ENV{TRACE_INC} )
find_library( IFBEAM               NAMES ifbeam            PATHS $ENV{IFBEAM_DIR} )

link_libraries( ${LIB_NAME} -L$ENV{BOOST_LIB} -lboost_system ${ROOTLIB} )

cet_make_library(LIBRARIES Boost::system
	LIBRARY_NAME sbn_BNBSpillInfoRetriever_MWRData
	SOURCE MWRData.cpp
)


simple_plugin(BNBRetriever module
	art::Persistency_Common
	art::Utilities canvas::canvas
	cetlib::cetlib cetlib_except::cetlib_except
	ROOT::X3d
        Boost::system
	messagefacility::MF_MessageLogger
	ifbeam::ifbeam
        ifdh_art::IFBeam_service
	sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_Common
  	sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays
	sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_ICARUS
        artdaq_core::artdaq-core_Utilities
	sbn_BNBSpillInfoRetriever_MWRData
        sbnobj::Common_POTAccounting
        larcorealg::CoreUtils
)

install_headers()
install_fhicl()
install_source()

