include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
include_directories( ${ROOT_INCLUDE_DIRS} )
include_directories( $ENV{NUTOOLS_INC} )
include_directories( $ENV{LARCOREOBJ_INC} )
include_directories( $ENV{LARCOREALG_INC} )
include_directories( $ENV{LARDATAALG_INC} )
include_directories( $ENV{LARDATAOBJ_INC} )
include_directories( $ENV{NUSIMDATA_INC} )

art_make(
  LIBRARY_NAME sbncode_DataTypes
  LIB_LIBRARIES
    ${ROOT_BASIC_LIB_LIST}
    ${ROOT_EG}
)

install_headers()
install_fhicl()
install_source()

