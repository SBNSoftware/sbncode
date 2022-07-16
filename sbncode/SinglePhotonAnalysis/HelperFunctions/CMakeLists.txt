art_make_library( LIBRARY_NAME HelperFunctions
                  SOURCE  helper_PandoraPFParticles.cxx
                  LIBRARIES 
                    ${ART_FRAMEWORK_CORE}
                    ${ART_FRAMEWORK_SERVICES_REGISTRY}
                    ${ART_FRAMEWORK_SERVICES}
                    ${ART_FRAMEWORK_PRINCIPAL}
                    art_Persistency_Common
                    art_Utilities canvas
                    ${ROOT_BASIC_LIB_LIST}
                    ${MF_MESSAGELOGGER}
                    ${MF_UTILITIES}
                    ${FHICLCPP}
                    cetlib cetlib_except
                    larsim_MCCheater_BackTrackerService_service
                    larsim_MCCheater_ParticleInventoryService_service
                    larcorealg_Geometry
                    )
