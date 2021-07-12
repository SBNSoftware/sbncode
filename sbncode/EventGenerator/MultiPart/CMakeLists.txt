# Add the subfolder containing the deconvolution tools
art_make( 
          MODULE_LIBRARIES larcorealg_Geometry
            larcore_Geometry_Geometry_service
            lardata_Utilities
            larevt_Filters
            lardataobj_RawData
            larevt_CalibrationDBI_IOVData
            larevt_CalibrationDBI_Providers
            lardataobj_RecoBase
            lardata_ArtDataHelper
            larsim_Simulation nug4_ParticleNavigation 
            lardataobj_Simulation
            nusimdata_SimulationBase
            lardata_Utilities
            lardataobj_MCBase
            larcoreobj_SummaryData
            nusimdata_SimulationBase
            nurandom_RandomUtils_NuRandomService_service
            lardataobj_MCBase
            ${CLHEP}
            ${ART_FRAMEWORK_CORE}
            ${ART_FRAMEWORK_PRINCIPAL}
            ${ART_FRAMEWORK_SERVICES_REGISTRY}
            ${ART_ROOT_IO_TFILE_SUPPORT} ${ROOT_CORE}
            ${ART_ROOT_IO_TFILESERVICE_SERVICE}
            art_Persistency_Common
            art_Persistency_Provenance
            art_Utilities
            canvas
            ${MF_MESSAGELOGGER}
            ${MF_UTILITIES}
            ${FHICLCPP}
            ${CETLIB}
            ${ROOT_FFTW}
            ${ROOT_BASIC_LIB_LIST}
	    ${ROOT_EG}
        )

install_headers()
install_fhicl()
install_source()

