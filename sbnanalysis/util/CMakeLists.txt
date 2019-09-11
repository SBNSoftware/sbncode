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
include_directories( $ENV{UBOONECODE_INC} )
include_directories($ENV{LARCOREALG_INC})
include_directories($ENV{CLHEP_INC})

link_directories(
  $ENV{CANVAS_LIB}
  $ENV{CETLIB_LIB}
  $ENV{GALLERY_LIB}
  $ENV{NUSIMDATA_LIB}
  $ENV{LARCOREOBJ_LIB}
  $ENV{CLHEP_LIB}
  $ENV{LARCOREALG_LIB}
  $ENV{LARDATAOBJ_LIB}
  ${CMAKE_INSTALL_DIR}/lib
)

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(sbnanalysis_util_Interaction Interaction.cxx)
target_link_libraries(sbnanalysis_util_Interaction ${ROOT_LIBRARIES})
install(TARGETS sbnanalysis_util_Interaction DESTINATION lib)

add_library(sbnanalysis_util_FakeDataSmearing FakeDataSmearing.cxx)
target_link_libraries(sbnanalysis_util_FakeDataSmearing ${ROOT_LIBRARIES})
install(TARGETS sbnanalysis_util_FakeDataSmearing DESTINATION lib)

