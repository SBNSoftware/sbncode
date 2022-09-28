add_subdirectory(TrackHit)
add_subdirectory(TrackSplit)
add_subdirectory(VertexStub)
add_subdirectory(CalorimetryAnalysis)
add_subdirectory(NCReco)

link_directories(
 $ENV{ROOTSYS}/lib
)

simple_plugin(NuVertexChargeTree module
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

	sbnobj_Common_Reco
	lardataobj_AnalysisBase
        
        art_root_io_tfile_support
        art_root_io
        art_root_io_dict
        caf_RecoUtils
        larsim_MCCheater_BackTrackerService_service
        nusimdata_SimulationBase
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
	larpandora_LArPandoraInterface
        larsim_MCCheater_ParticleInventoryService_service
	sbn_LArReco	    
	sbncode_TPCReco_VertexStub
)
simple_plugin ( EventSelect module
                IFDH_SERVICE
               ${ART_FRAMEWORK_CORE}
               ${ART_FRAMEWORK_PRINCIPAL}
               ${ART_FRAMEWORK_SERVICES_REGISTRY}
               ${ART_FRAMEWORK_SERVICES_BASIC}
               ${ART_FRAMEWORK_SERVICES_OPTIONAL}
               art_Utilities canvas
               ${MF_MESSAGELOGGER}
               ${MF_UTILITIES}
               ${FHICLCPP}
               cetlib cetlib_except
	${ROOT_EVE_LIB_LIST}
	${ROOT_X3d}
	${ROOT_BASIC_LIB_LIST}
)


install_headers()
install_fhicl()
install_source()

