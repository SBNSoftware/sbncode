# for classes_def.xml!!
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( $ENV{NUTOOLS_INC} )

art_make( LIBRARY_NAME sbn_TPCReco_Products
          LIB_LIBRARIES
                        ${ROOT_BASIC_LIB_LIST}
                        ${ROOT_EG}
          DICT_LIBRARIES
                        ${ROOT_BASIC_LIB_LIST}
                        ${ROOT_EG}
        )

install_headers()
install_fhicl()
install_source()
