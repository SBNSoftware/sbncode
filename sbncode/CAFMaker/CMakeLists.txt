# for classes_def.xml
include_directories ( ${CMAKE_CURRENT_SOURCE_DIR} )
set( PACKAGE CAFMaker )
FILE( GLOB src_files *.cxx )

add_subdirectory(RecoUtils)

art_make_library( LIBRARY_NAME ${PACKAGE}
                  SOURCE ${src_files}
                  LIBRARIES 
                            StandardRecord
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
	art_Framework_Services_System_TriggerNamesService_service
                )

SET(ENV{IFHDC_LIB}   "$ENV{IFDHC_FQ_DIR}/lib")
cet_find_library( IFDH NAMES ifdh PATHS ENV IFDHC_LIB NO_DEFAULT_PATH )
include_directories( $ENV{IFDHC_FQ_DIR}/inc )
include_directories( $ENV{IFDH_ART_FQ_DIR}/inc )

cet_find_library( IFDH_SERVICE NAMES IFDH_service PATHS ENV IFDH_ART_LIB )

simple_plugin ( CAFMaker module
                CAFMaker
               StandardRecord
 		            ${IFDH_SERVICE}
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
               ${ROOT_BASIC_LIB_LIST}
		art_root_io_RootDB
		hep_concurrency
		nurandom_RandomUtils_NuRandomService_service
               BASENAME_ONLY
            )

build_dictionary( CAFMaker
                  DICTIONARY_LIBRARIES 
				StandardRecord
				StandardRecord_dict
				${ART_FRAMEWORK_CORE}
				${ART_FRAMEWORK_PRINCIPAL}
				art_Persistency_Common
				art_Persistency_Provenance
				art_Utilities canvas
				cetlib cetlib_except
				${MF_MESSAGELOGGER}
				${MF_UTILITIES}
)


install_headers()
install_fhicl()
install_source()
