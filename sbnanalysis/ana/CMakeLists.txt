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

link_directories(
  $ENV{CANVAS_LIB}
  $ENV{CETLIB_LIB}
  $ENV{GALLERY_LIB}
  $ENV{NUSIMDATA_LIB}
  $ENV{LARCOREOBJ_LIB}
  $ENV{LARDATAOBJ_LIB}
)

if (DEFINED ENV{SBNANA_STATIC_PROCESSOR})
  add_definitions(-DCMAKE_PROCESSOR_INCLUDE=\"$ENV{SBNANA_STATIC_PROCESSOR}.h\")

  Message(STATUS ${CMAKE_INSTALL_PREFIX}/lib/libsbnanalysis_ExampleAnalysis_ExampleSelection.so)
  include_directories( ${CMAKE_INSTALL_PREFIX}/lib )
  link_directories(
    $ENV{JSONCPP_LIB}
    ${CMAKE_INSTALL_PREFIX}/lib
  )
  add_executable(static static.cxx)
  target_link_libraries(static sbnanalysis_Processor jsoncpp sbnanalysis_ExampleAnalysis_ExampleSelection)
  install(TARGETS static DESTINATION bin)
endif()


add_subdirectory(ExampleAnalysis)



