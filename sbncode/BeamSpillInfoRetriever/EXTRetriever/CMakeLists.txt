
include_directories( $ENV{IFBEAM_FQ_DIR}/include )
include_directories( $ENV{LIBWDA_FQ_DIR}/include )
include_directories( $ENV{IFDHC_FQ_DIR}/inc )
include_directories( $ENV{IFDH_ART_INC} )

simple_plugin(EXTRetriever module
	${ART_FRAMEWORK_CORE}
	${ART_FRAMEWORK_SERVICES_REGISTRY}
	${ART_FRAMEWORK_SERVICES}
	${ART_FRAMEWORK_PRINCIPAL}
	art_Persistency_Common
	art_Utilities canvas
	${FHICLCPP}
	cetlib cetlib_except
	MF_MessageLogger
       	sbnobj_Common_POTAccounting
	sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_Common
  	sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays
	sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_ICARUS	
        artdaq_core::artdaq-core_Utilities
)

install_headers()
install_fhicl()
install_source()

