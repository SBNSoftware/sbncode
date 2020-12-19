include_directories( $ENV{IFDHC_LIB}/../inc )
link_directories( $ENV{ROOTSYS}/lib )


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
)

simple_plugin( NuMiKaonGen tool
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
)

add_subdirectory(Higgs)
add_subdirectory(HNL)

install_headers()
install_fhicl()
install_source()
