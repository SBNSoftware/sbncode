link_directories(
 $ENV{ROOTSYS}/lib
)

add_subdirectory(LArReco)
add_subdirectory(Products)

simple_plugin(MCSFitAllPID module
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
	${ROOT_X3d}
	${ROOT_BASIC_LIB_LIST}
	MF_MessageLogger

	sbn_LArReco	    
	sbncode_LArRecoProducer_Products_dict
	lardataobj_RecoBase_AssnsDicts_dict
	lardataobj_RecoBase_dict	
)

simple_plugin(RangePAllPID module
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
	${ROOT_X3d}
	${ROOT_BASIC_LIB_LIST}
	MF_MessageLogger

	sbn_LArReco	    
	sbncode_LArRecoProducer_Products_dict
)

install_headers()
install_fhicl()
install_source()

