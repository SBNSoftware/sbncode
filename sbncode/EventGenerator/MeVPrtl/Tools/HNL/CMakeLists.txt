art_make_library( LIBRARY_NAME sbncode_EventGenerator_MeVPrtl_Dalitz
          SOURCE HNLDecayDalitz.cpp
          LIBRARIES
                        ${ROOT_BASIC_LIB_LIST}
                        ${ROOT_EVE_LIB_LIST}
                        ${ROOT_X3d}
			sbncode_EventGenerator_MeVPrtl_Constants
        )

simple_plugin( Kaon2HNLFlux tool
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
			EGPythia6
			${CLHEP}
			sbnobj_Common_EventGen_MeVPrtl
			nugen_EventGeneratorBase_GENIE
			nurandom_RandomUtils_NuRandomService_service
			gsl
			sbncode_EventGenerator_MeVPrtl_Constants
)

simple_plugin( HNLMakeDecay tool
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
			${CLHEP}
			nurandom_RandomUtils_NuRandomService_service
			nusimdata_SimulationBase
			gsl
			sbncode_EventGenerator_MeVPrtl_Dalitz
			sbncode_EventGenerator_MeVPrtl_Constants
)


install_headers()
install_fhicl()
install_source()
