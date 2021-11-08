
art_make_library( LIBRARY_NAME sbncode_EventGenerator_MeVPrtl_Constants
          SOURCE Constants.cpp
          LIBRARIES
                        ${ROOT_BASIC_LIB_LIST}
                        ${ROOT_EVE_LIB_LIST}
                            ${FHICLCPP}
                            cetlib cetlib_except
        )


simple_plugin( EmptyKaonGen tool
                            ${ART_FRAMEWORK_CORE}
                            ${ART_FRAMEWORK_SERVICES_REGISTRY}
                            ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                            ${ART_FRAMEWORK_SERVICES}
                            ${ART_FRAMEWORK_PRINCIPAL}
                            art_Persistency_Common
                            art_Utilities canvas
                            ${MF_MESSAGELOGGER}
                            ${MF_UTILITIES}
                            ${FHICLCPP}
                            cetlib cetlib_except
                            ${ROOT_EVE_LIB_LIST}
                            ${ROOT_BASIC_LIB_LIST}
			${CLHEP}
			nurandom_RandomUtils_NuRandomService_service
			nusimdata_SimulationBase
			sbncode_EventGenerator_MeVPrtl_Constants
)

simple_plugin( NuMiKaonGen tool
                            IFDH_SERVICE
                            ${ART_FRAMEWORK_CORE}
                            ${ART_FRAMEWORK_SERVICES_REGISTRY}
                            ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                            ${ART_FRAMEWORK_SERVICES}
                            ${ART_FRAMEWORK_PRINCIPAL}
                            art_Persistency_Common
                            art_Utilities canvas
                            ${MF_MESSAGELOGGER}
                            ${MF_UTILITIES}
                            ${FHICLCPP}
                            cetlib cetlib_except
                            ${ROOT_EVE_LIB_LIST}
                            ${ROOT_BASIC_LIB_LIST}
			${CLHEP}
			ifdh
			dk2nuTree
			nurandom_RandomUtils_NuRandomService_service
			nusimdata_SimulationBase
			sbncode_EventGenerator_MeVPrtl_Constants
)

simple_plugin( WeightedRayTraceBox tool
                            ${ART_FRAMEWORK_CORE}
                            ${ART_FRAMEWORK_SERVICES_REGISTRY}
                            ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                            ${ART_FRAMEWORK_SERVICES}
                            ${ART_FRAMEWORK_PRINCIPAL}
                            art_Persistency_Common
                            art_Utilities canvas
                            ${MF_MESSAGELOGGER}
                            ${MF_UTILITIES}
                            ${FHICLCPP}
                            cetlib cetlib_except
                            ${ROOT_EVE_LIB_LIST}
                            ${ROOT_BASIC_LIB_LIST}
			${CLHEP}
			larcorealg_Geometry
			nurandom_RandomUtils_NuRandomService_service
			dk2nuGenie
			dk2nuTree
			sbncode_EventGenerator_MeVPrtl_Constants
)

simple_plugin( MixedWeightRayTraceBox tool
                            ${ART_FRAMEWORK_CORE}
                            ${ART_FRAMEWORK_SERVICES_REGISTRY}
                            ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                            ${ART_FRAMEWORK_SERVICES}
                            ${ART_FRAMEWORK_PRINCIPAL}
                            art_Persistency_Common
                            art_Utilities canvas
                            ${MF_MESSAGELOGGER}
                            ${MF_UTILITIES}
                            ${FHICLCPP}
                            cetlib cetlib_except
                            ${ROOT_EVE_LIB_LIST}
                            ${ROOT_BASIC_LIB_LIST}
			${CLHEP}
			larcorealg_Geometry
			nurandom_RandomUtils_NuRandomService_service
			dk2nuTree
			sbncode_EventGenerator_MeVPrtl_Constants
)




simple_plugin( ReThrowRayTraceBox tool
                            ${ART_FRAMEWORK_CORE}
                            ${ART_FRAMEWORK_SERVICES_REGISTRY}
                            ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                            ${ART_FRAMEWORK_SERVICES}
                            ${ART_FRAMEWORK_PRINCIPAL}
                            art_Persistency_Common
                            art_Utilities canvas
                            ${MF_MESSAGELOGGER}
                            ${MF_UTILITIES}
                            ${FHICLCPP}
                            cetlib cetlib_except
                            ${ROOT_EVE_LIB_LIST}
                            ${ROOT_BASIC_LIB_LIST}
			${CLHEP}
			larcorealg_Geometry
			nurandom_RandomUtils_NuRandomService_service
			sbncode_EventGenerator_MeVPrtl_Constants
)

simple_plugin( RayTraceBox tool
                            ${ART_FRAMEWORK_CORE}
                            ${ART_FRAMEWORK_SERVICES_REGISTRY}
                            ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                            ${ART_FRAMEWORK_SERVICES}
                            ${ART_FRAMEWORK_PRINCIPAL}
                            art_Persistency_Common
                            art_Utilities canvas
                            ${MF_MESSAGELOGGER}
                            ${MF_UTILITIES}
                            ${FHICLCPP}
                            cetlib cetlib_except
                            ${ROOT_EVE_LIB_LIST}
                            ${ROOT_BASIC_LIB_LIST}
			${CLHEP}
			larcorealg_Geometry
			nurandom_RandomUtils_NuRandomService_service
			sbncode_EventGenerator_MeVPrtl_Constants
)

add_subdirectory(Higgs)
add_subdirectory(HNL)

install_headers()
install_fhicl()
install_source()
