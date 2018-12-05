include_directories( ${ROOT_INCLUDE_DIRS} )
include_directories( $ENV{BOOST_INC})
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
include_directories( $ENV{UBCORE_INC} )
include_directories( $ENV{CLHEP_INC} )
include_directories( $ENV{LARDATA_INC} )
include_directories( $ENV{LARRECO_INC} )
include_directories( $ENV{LARCOREALG_INC} )
include_directories( $ENV{LARDATAALG_INC} )
include_directories( $ENV{ICARUSCODE_INC} )
include_directories( $ENV{SBNDCODE_INC} )
include_directories( $ENV{ICARUSCODE_INC} )
include_directories( $ENV{MESSAGEFACILITY_INC} )

link_directories(
  $ENV{MESSAGEFACILITY_LIB}
  $ENV{FHICLCPP_LIB}
  $ENV{BOOST_LIB}
  $ENV{CANVAS_LIB}
  $ENV{CANVAS_ROOT_IO_LIB}
  $ENV{CETLIB_EXCEPT_LIB}
  $ENV{CETLIB_LIB}
  $ENV{GALLERY_LIB}
  $ENV{NUSIMDATA_LIB}
  $ENV{LARCOREOBJ_LIB}
  $ENV{LARDATAOBJ_LIB}
  $ENV{UBCORE_LIB}
  $ENV{LARDATAALG_LIB}
  $ENV{LARSIM_LIB}
  $ENV{LARCOREALG_LIB}
  $ENV{LARCOREALG_LIB}
  $ENV{LARDATA_LIB}
  $ENV{LARRECO_LIB}
  $ENV{ICARUSCODE_LIB}
  $ENV{SBNDCODE_LIB}
  $ENV{LARSIM_LIB}
  $ENV{ICARUSCODE_LIB}
  $ENV{SBNDCODE_LIB}
  $ENV{JSONCPP_LIB}
  $ENV{UBCORE_LIB}
  ${CMAKE_INSTALL_DIR}/lib
)

add_subdirectory(ExampleAnalysis)
add_subdirectory(SelectionAnalysis)
add_subdirectory(SBNOsc)

