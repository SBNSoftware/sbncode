simple_plugin( CalorimetryAnalysis module
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
	art_root_io_tfile_support
	art_root_io
	art_root_io_dict
        caf_RecoUtils
        lardataobj_AnalysisBase
        lardataobj_RecoBase
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
        larcorealg_GeoAlgo
        larsim_MCCheater_BackTrackerService_service
        nusimdata_SimulationBase
        larsim_MCCheater_ParticleInventoryService_service
        sbncode_LArRecoProducer_Products_dict
	larcorealg_GeoAlgo
)
install_headers()
install_fhicl()
install_source()

