link_directories(
 $ENV{ROOTSYS}/lib
)

add_subdirectory(LArReco)
add_subdirectory(Products)

art_make_library(
	LIBRARY_NAME sbn_LArReco_PCA
	SOURCE PCA.cc
	LIBRARIES
	sbncode_LArRecoProducer_Products_dict
	lardataobj_RecoBase_AssnsDicts_dict
	lardataobj_RecoBase_dict	
	sbndcode_OpDetSim
)

simple_plugin(PMTFlashTriggerFilter module
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
	sbncode_LArRecoProducer_Products_dict
	lardataobj_RecoBase_AssnsDicts_dict
	lardataobj_RecoBase_dict	
        larcorealg_Geometry
	sbn_LArReco	    
)


simple_plugin(PMTFlashTriggerMaker module
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
	sbndcode_OpDetSim
)

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
	lardataobj_RecoBase
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

simple_plugin(TrackAreaHit module
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
        lardataobj_RecoBase_AssnsDicts_dict
        lardataobj_RecoBase     
        lardataobj_RecoBase_dict        

        larcorealg_Geometry
        larcore_Geometry_Geometry_service
	lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
	lardata_DetectorInfoServices_DetectorPropertiesServiceStandard_service
	lardataalg_DetectorInfo
)

simple_plugin(TrackHitFilter module
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
        lardataobj_RecoBase_AssnsDicts_dict
        lardataobj_RecoBase     
        lardataobj_RecoBase_dict        

        larcorealg_Geometry
        larcore_Geometry_Geometry_service
)

simple_plugin(TrackHitDumper module
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
	lardataobj_RecoBase_AssnsDicts_dict
	lardataobj_RecoBase	
	lardataobj_RecoBase_dict	

	larcorealg_Geometry
	larcore_Geometry_Geometry_service
)

simple_plugin(ShowerSelectionVars module
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

simple_plugin(TrackSplitter module
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
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  lardata_DetectorInfoServices_DetectorClocksServiceStandard_service
  lardata_DetectorInfoServices_DetectorPropertiesServiceStandard_service

  sbn_LArReco_PCA
  sbn_LArReco
  sbncode_LArRecoProducer_Products_dict
)

simple_plugin(MergedTrackIdentifier module
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
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  larpandora_LArPandoraInterface

  sbn_LArReco
  sbncode_LArRecoProducer_Products_dict
)



simple_plugin(CRTHitConverter module
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

	sbncode_LArRecoProducer_Products_dict
	lardataobj_RecoBase_AssnsDicts_dict
	lardataobj_RecoBase_dict	
	sbndcode_CRT_CRTProducts_dict
	icaruscode_CRTProducts	
)

simple_plugin(CRTHitMatch module
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

	lardataobj_RecoBase
	sbncode_LArRecoProducer_Products_dict
	sbndcode_CRT_CRTProducts_dict
	sbndcode_CRTUtils
	lardataobj_AnalysisBase
)

simple_plugin(VertexChargeVacuum module
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

	sbncode_LArRecoProducer_Products_dict
	lardataobj_AnalysisBase
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
	sbn_LArReco_PCA
)

simple_plugin(PCAnglePlaneMaker module
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

	sbncode_LArRecoProducer_Products_dict
	lardataobj_AnalysisBase
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
	sbn_LArReco_PCA
)

simple_plugin(PCAngleKinkFinder module
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
	Minuit

	sbncode_LArRecoProducer_Products_dict
	lardataobj_AnalysisBase
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
)

simple_plugin(PCAngleKinkTree module
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

	sbncode_LArRecoProducer_Products_dict
	lardataobj_AnalysisBase
        
        art_root_io_tfile_support
        art_root_io
        art_root_io_dict
        caf_RecoUtils
        larsim_MCCheater_BackTrackerService_service
        nusimdata_SimulationBase
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
	sbn_LArReco_PCA
)

install_headers()
install_fhicl()
install_source()

