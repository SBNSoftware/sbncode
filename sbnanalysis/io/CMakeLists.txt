include_directories( ${CMAKE_SOURCE_DIR} )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

include_directories( ${ROOT_INCLUDE_DIRS} )
include_directories( $ENV{NUTOOLS_INC} )
include_directories( $ENV{FHICLCPP_INC} )
include_directories( $ENV{CANVAS_INC} )
include_directories( $ENV{CANVAS_ROOT_IO_INC} )
include_directories( $ENV{CETLIB_INC} )
include_directories( $ENV{CETLIB_EXCEPT_INC} )
include_directories( $ENV{GALLERY_INC} )
include_directories( $ENV{LARCOREOBJ_INC} )
include_directories( $ENV{LARDATAOBJ_INC} )
include_directories( $ENV{NUSIMDATA_INC} )

link_directories(
  $ENV{CANVAS_LIB}
  $ENV{CETLIB_LIB}
  $ENV{GALLERY_LIB}
  $ENV{NUSIMDATA_LIB}
  $ENV{LARCOREOBJ_LIB}
  $ENV{LARDATAOBJ_LIB}
)

ROOT_GENERATE_DICTIONARY(sbn_eventdict Event.hh LINKDEF linkdef.h)

add_library(sbnanalysis_Processor SHARED ProcessorBase.cxx SelectionBase.cxx sbn_eventdict.cxx)

target_link_libraries(sbnanalysis_Processor ${ROOT_LIBRARIES})

file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/lib)
add_custom_command(TARGET sbnanalysis_Processor POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
  libsbn_eventdict_rdict.pcm ${CMAKE_INSTALL_PREFIX}/lib)

install(TARGETS sbnanalysis_Processor DESTINATION lib)

