# for classes_def.xml!!
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )
# include_directories( $ENV{NUTOOLS_INC} )
get_target_property(DICT_INC_DIR lardataobj::RecoBase INTERFACE_INCLUDE_DIRECTORIES)
include_directories( ${DICT_INC_DIR} )
get_target_property(DICT_INC_DIR larcorealg::Geometry INTERFACE_INCLUDE_DIRECTORIES)
include_directories( ${DICT_INC_DIR} )
get_target_property(DICT_INC_DIR larcoreobj::SimpleTypesAndConstants INTERFACE_INCLUDE_DIRECTORIES)
include_directories( ${DICT_INC_DIR} )

art_make(
    LIB_LIBRARIES
        lardataobj::RecoBase
        ROOT::Core
        ROOT::Eg
)

install_headers()
install_fhicl()
install_source()
