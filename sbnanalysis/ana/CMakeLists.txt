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
include_directories( $ENV{JSONCPP_INC} )
# include stuff from uboonecode
include_directories( $ENV{UBOONECODE_INC} )

# determine location of uboonecode
if (DEFINED ENV{UBOONECODE_LIB})
  set(UBOONE_PATH $ENV{UBOONECODE_LIB})
else()
  set(UBOONE_PATH $ENV{MRB_BUILDDIR}/uboonecode/lib/)
endif()

link_directories(
  $ENV{CANVAS_LIB}
  $ENV{CETLIB_LIB}
  $ENV{GALLERY_LIB}
  $ENV{NUSIMDATA_LIB}
  $ENV{LARCOREOBJ_LIB}
  $ENV{LARDATAOBJ_LIB}
  ${UBOONE_PATH}
  ${CMAKE_INSTALL_DIR}/lib
)

add_subdirectory(ExampleAnalysis)
add_subdirectory(SBNOsc)

