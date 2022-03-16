
include_directories( $ENV{IFBEAM_FQ_DIR}/include )
include_directories( $ENV{LIBWDA_FQ_DIR}/include )
include_directories( $ENV{IFDHC_FQ_DIR}/inc )
include_directories( $ENV{IFDH_ART_INC} )

simple_plugin(EXTRetriever module
	art::Framework_Core
	art::Framework_Services_Registry
	${ART_FRAMEWORK_SERVICES}
	art::Framework_Principal
	art::Persistency_Common
	art::Utilities canvas::canvas
	fhiclcpp::fhiclcpp
	cetlib::cetlib cetlib_except
	messagefacility::MF_MessageLogger
       	sbnobj_Common_POTAccounting
	sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_Common
  	sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays
	sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_ICARUS	
        artdaq_core::artdaq-core_Utilities
)

install_headers()
install_fhicl()
install_source()

