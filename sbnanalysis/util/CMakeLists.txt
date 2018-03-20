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

link_directories(
  $ENV{CANVAS_LIB}
  $ENV{CETLIB_LIB}
  $ENV{GALLERY_LIB}
  $ENV{NUSIMDATA_LIB}
  $ENV{LARCOREOBJ_LIB}
  $ENV{LARDATAOBJ_LIB}
  ${CMAKE_INSTALL_DIR}/lib
)

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(sbnanalysis_util_Covariance Covariance.cxx)
target_link_libraries(sbnanalysis_util_Covariance ${ROOT_LIBRARIES})

add_executable(sbn-cov sbn_cov.cxx)
target_link_libraries(sbn-cov sbnanalysis_util_Covariance sbnanalysis_Event)

install(TARGETS sbnanalysis_util_Covariance DESTINATION lib)
install(TARGETS sbn-cov DESTINATION bin)

