find_ups_product(ifdhc)
find_ups_product(dk2nudata)
find_ups_product(dk2nugenie)

link_directories( $ENV{DK2NUDATA_LIB} )
link_directories( $ENV{DK2NUGENIE_LIB} )

add_subdirectory(Tools)
add_subdirectory(config)

simple_plugin( MeVPrtlGen module
              ${ART_FRAMEWORK_CORE}
              ${ART_FRAMEWORK_SERVICES_REGISTRY}
              ${ART_ROOT_IO_TFILESERVICE_SERVICE}
              ${ART_FRAMEWORK_SERVICES}
				      ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
				      ${ART_ROOT_IO_TFILE_SUPPORT}
              ${ART_FRAMEWORK_PRINCIPAL}
              art_Persistency_Common
              art_Utilities canvas
              ${MF_MESSAGELOGGER}
              ${MF_UTILITIES}
              ${FHICLCPP}
              cetlib cetlib_except
              ${ROOT_EVE_LIB_LIST}
              ${ROOT_X3d}
              ${ROOT_CORE}
              ${ROOT_BASIC_LIB_LIST}
        			${ROOT_XMLIO}
        			${CLHEP}
		nusimdata_SimulationBase
              larcorealg_Geometry
        			nurandom_RandomUtils_NuRandomService_service
        			larcoreobj_SummaryData
        			sbnobj_Common_EventGen_MeVPrtl
        			sbncode_EventGenerator_MeVPrtl_Constants
)

simple_plugin( MeVPrtlTestRayTrace module
                            ${ART_FRAMEWORK_CORE}
                            ${ART_FRAMEWORK_SERVICES_REGISTRY}
                            ${ART_ROOT_IO_TFILESERVICE_SERVICE}
                            ${ART_FRAMEWORK_SERVICES}
				${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
				${ART_ROOT_IO_TFILE_SUPPORT}
                            ${ART_FRAMEWORK_PRINCIPAL}
                            art_Persistency_Common
                            art_Utilities canvas
                            ${MF_MESSAGELOGGER}
                            ${MF_UTILITIES}
                            ${FHICLCPP}
                            cetlib cetlib_except
                            ${ROOT_EVE_LIB_LIST}
                            ${ROOT_X3d}
                            ${ROOT_CORE}
                            ${ROOT_BASIC_LIB_LIST}
			${ROOT_XMLIO}
			${CLHEP}
			nurandom_RandomUtils_NuRandomService_service
			larcoreobj_SummaryData
        		sbnobj_Common_EventGen_MeVPrtl
			sbncode_EventGenerator_MeVPrtl_Constants
)



install_headers()
install_fhicl()
install_source()
