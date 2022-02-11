link_libraries( -L$ENV{BOOST_LIB} -lboost_system -lboost_filesystem)

art_make( BASENAME_ONLY
          LIB_LIBRARIES
            larcorealg_Geometry
            larcore_Geometry_Geometry_service
            lardata_Utilities
            larevt_Filters
            larpandoraobj
            lardataobj_RawData
			lardataobj_RecoBase
            lardataobj_MCBase
			larcoreobj_SummaryData
			larsim_MCCheater_BackTrackerService_service
			larsim_MCCheater_ParticleInventoryService_service
			larpandora_LArPandoraInterface
            ${PANDORASDK}
            ${PANDORAMONITORING}
            LArPandoraContent 
			nusimdata_SimulationBase
			${Boost_SERIALIZATION_LIBRARY}
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        
                        ${FHICLCPP}
			cetlib cetlib_except
                        ${CLHEP}
                        ${ROOT_BASIC_LIB_LIST}
                        ${ROOT_GEOM} 
          MODULE_LIBRARIES
                        larcorealg_Geometry
                        larcore_Geometry_Geometry_service
                        lardata_Utilities
                        larevt_Filters
                        lardataobj_RawData
			lardataobj_RecoBase
            lardataobj_MCBase
            larcoreobj_SummaryData	
			larsim_MCCheater_BackTrackerService_service
			larpandora_LArPandoraInterface
			nusimdata_SimulationBase
			${Boost_SERIALIZATION_LIBRARY}
                        ${ART_FRAMEWORK_CORE}
                        ${ART_FRAMEWORK_PRINCIPAL}
                        ${ART_FRAMEWORK_SERVICES_REGISTRY}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_RANDOMNUMBERGENERATOR_SERVICE}
                        ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
                        art_Persistency_Common canvas
                        art_Persistency_Provenance canvas
                        art_Utilities canvas
                        ${MF_MESSAGELOGGER}
                        
                        ${FHICLCPP}
			cetlib cetlib_except
                        ${CLHEP}
                        ${ROOT_BASIC_LIB_LIST}
                        ${ROOT_GEOM} )
                      
#add_subdirectory(SEAview)
add_subdirectory(job)
install_headers()
install_fhicl()
install_source()
