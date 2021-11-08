
art_make_library(
	LIBRARY_NAME sbn_BNBSpillInfoRetriever_MWRData
	SOURCE MWRData.cpp
)


simple_plugin(BNBRetriever module
	${ART_FRAMEWORK_CORE}
	${ART_FRAMEWORK_SERVICES_REGISTRY}
	${ART_ROOT_IO_TFILESERVICE_SERVICE}                ${ART_FRAMEWORK_SERVICES}
	${ART_FRAMEWORK_PRINCIPAL}
	art_Persistency_Common
	art_Utilities canvas
	${MF_MESSAGELOGGER}
	${MF_UTILITIES}
	${FHICLCPP}
	cetlib cetlib_except
	${ROOT_EVE_LIB_LIST}
	${ROOT_BASIC_LIB_LIST}
        ${Boost_PROGRAM_OPTIONS_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} ${Boost_DATE_TIME_LIBRARY} ${Boost_SERIALIZATION_LIBRARY}
	MF_MessageLogger
       	sbnobj_Common_POTAccounting
	#${IFBEAMSERVICE}
	#${IFBEAM}
	ifbeam
	ifdh_art::IFBeam_service
	sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_Common
  	sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays
	sbndaq_artdaq_core::sbndaq-artdaq-core_Overlays_ICARUS	
        artdaq_core::artdaq-core_Utilities
	sbn_BNBSpillInfoRetriever_MWRData	
)

install_headers()
install_fhicl()
install_source()

