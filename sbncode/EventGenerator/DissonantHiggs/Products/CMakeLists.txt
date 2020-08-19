# for classes_def.xml!!
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( $ENV{NUTOOLS_INC} )

FILE( GLOB src_files *.cxx )

art_make_library( LIBRARY_NAME sbncode_EventGenerator_DissonantHiggs_Products
          SOURCE ${src_files}
          LIBRARIES
                        ${ROOT_BASIC_LIB_LIST}
			nusimdata_SimulationBase
        )

build_dictionary(sbncode_EventGenerator_DissonantHiggs_Products 
                  DICTIONARY_LIBRARIES sbncode_EventGenerator_DissonantHiggs_Products
                  ${ART_FRAMEWORK_CORE}
                                  ${ART_FRAMEWORK_PRINCIPAL}
                                  art_Persistency_Common
                                  art_Utilities canvas
                                  cetlib cetlib_except
                                  ${MF_MESSAGELOGGER}
                                  ${MF_UTILITIES}
                 )


install_headers()
install_fhicl()
install_source()
