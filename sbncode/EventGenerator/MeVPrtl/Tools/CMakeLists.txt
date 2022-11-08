link_directories( $ENV{ROOTSYS}/lib )

art_make_library( LIBRARY_NAME sbncode_EventGenerator_MeVPrtl_Constants
          SOURCE Constants.cpp
          LIBRARIES
                        ${ROOT_BASIC_LIB_LIST}
                        ${ROOT_EVE_LIB_LIST}
                        ${ROOT_X3d}
                            ${FHICLCPP}
                            cetlib cetlib_except
        )


simple_plugin( EmptyKaonGen tool
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
			sbncode_EventGenerator_MeVPrtl_Constants
)

cet_find_library( IFDH NAMES ifdh PATHS ENV IFDHC_LIB NO_DEFAULT_PATH )
##include_directories( $ENV{IFDHC_FQ_DIR}/inc )
##include_directories( $ENV{IFDH_ART_FQ_DIR}/inc )

cet_find_library( IFDH_SERVICE NAMES IFDH_service PATHS ENV IFDH_ART_LIB )


simple_plugin( NuMiKaonGen tool
                            IFDH_SERVICE
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
			ifdh
			dk2nuTree
			nurandom_RandomUtils_NuRandomService_service
			nusimdata_SimulationBase
			sbncode_EventGenerator_MeVPrtl_Constants
)

simple_plugin( TxtFileGen tool
                            IFDH_SERVICE
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
			ifdh
			nurandom_RandomUtils_NuRandomService_service
			dk2nuTree
			nusimdata_SimulationBase
)

simple_plugin( WeightedRayTraceBox tool
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
			larcorealg_Geometry
			nurandom_RandomUtils_NuRandomService_service
			dk2nuGenie
			dk2nuTree
			sbncode_EventGenerator_MeVPrtl_Constants
)

simple_plugin( MixedWeightRayTraceBox tool
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
			larcorealg_Geometry
			nurandom_RandomUtils_NuRandomService_service
			dk2nuTree
			sbncode_EventGenerator_MeVPrtl_Constants
)




simple_plugin( ReThrowRayTraceBox tool
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
			larcorealg_Geometry
			nurandom_RandomUtils_NuRandomService_service
			sbncode_EventGenerator_MeVPrtl_Constants
)

simple_plugin( RayTraceBox tool
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
			larcorealg_Geometry
			nurandom_RandomUtils_NuRandomService_service
			sbncode_EventGenerator_MeVPrtl_Constants
)

add_subdirectory(Higgs)
add_subdirectory(HNL)
add_subdirectory(ALP)

install_headers()
install_fhicl()
install_source()
