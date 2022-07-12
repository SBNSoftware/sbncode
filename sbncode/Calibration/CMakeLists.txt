link_directories(
 $ENV{ROOTSYS}/lib
)

simple_plugin( TrackCaloSkimmer module
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
        lardataobj_AnalysisBase
        lardataobj_RecoBase
        larcorealg_Geometry
        larcore_Geometry_Geometry_service
        larcorealg_GeoAlgo
	sbnobj_Common_Reco
	larcorealg_GeoAlgo
	lardata_DetectorInfoServices_DetectorPropertiesServiceStandard_service
	lardataalg_DetectorInfo
	Minuit
	sbncode_CAFMaker
)

simple_plugin(TrackCaloSkimmerSelectStoppingTrack tool
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
        larcorealg_Geometry
	larcorealg_GeoAlgo
	lardata_DetectorInfoServices_DetectorPropertiesServiceStandard_service
	lardataalg_DetectorInfo
)

simple_plugin(TrackCaloSkimmerSelectAnode2CathodeTrack tool
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
)

simple_plugin ( DataSelect module
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


simple_plugin ( AnodeCrossT0 module
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
        lardataobj_AnalysisBase
        lardataobj_RecoBase
        larcorealg_Geometry
	larcorealg_GeoAlgo
	lardata_DetectorInfoServices_DetectorPropertiesServiceStandard_service
	lardataalg_DetectorInfo
)

add_subdirectory(fcl)

install_headers()
install_fhicl()
install_source()

