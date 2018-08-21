include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

include_directories(${ROOT_INCLUDE_DIRS})
include_directories($ENV{ART_INC})
include_directories($ENV{NUTOOLS_INC})
include_directories($ENV{FHICLCPP_INC})
include_directories($ENV{CANVAS_INC})
include_directories($ENV{CANVAS_ROOT_IO_INC})
include_directories($ENV{CETLIB_INC})
include_directories($ENV{CETLIB_EXCEPT_INC})
include_directories($ENV{GALLERY_INC})
include_directories($ENV{LARCOREOBJ_INC})
include_directories($ENV{LARDATAOBJ_INC})
include_directories($ENV{NUSIMDATA_INC})
include_directories($ENV{JSONCPP_INC})
include_directories($ENV{LARSIM_INC})

link_directories(
  $ENV{CANVAS_LIB}
  $ENV{CETLIB_LIB}
  $ENV{GALLERY_LIB}
  $ENV{NUSIMDATA_LIB}
  $ENV{LARCOREOBJ_LIB}
  $ENV{LARDATAOBJ_LIB}
  $ENV{JSONCPP_LIB}
  $ENV{LARSIM_LIB}
)

ROOT_GENERATE_DICTIONARY(sbn_eventdict Event.hh LINKDEF linkdef.h)

find_library(LARSIM_BASE_LIBRARY larsim_EventWeight_Base HINTS $ENV{LARSIM_LIB})
find_library(LARSIM_BASE_DICT larsim_EventWeight_Base_dict HINTS $ENV{LARSIM_LIB})

add_library(sbnanalysis_Event SHARED sbn_eventdict.cxx Event.cxx)
target_link_libraries(
  sbnanalysis_Event
  ${LARSIM_BASE_LIBRARY}
  ${LARSIM_BASE_DICT}
  ${ROOT_LIBRARIES}
)

add_library(sbnanalysis_Processor SHARED ProcessorBase.cxx SelectionBase.cxx ProcessorBlock.cxx Loader.cxx)
target_link_libraries(
  sbnanalysis_Processor
  sbnanalysis_Event
  canvas
  cetlib
  gallery
  nusimdata_SimulationBase
  larcoreobj_SummaryData
  lardataobj_RecoBase
  lardataobj_MCBase
  ${LARSIM_BASE_LIBRARY}
  ${LARSIM_BASE_DICT}
  ${ROOT_LIBRARIES}
)

file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/lib)
add_custom_command(TARGET sbnanalysis_Event POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
  libsbn_eventdict_rdict.pcm ${CMAKE_INSTALL_PREFIX}/lib)

add_executable(sbn sbn.cxx)
target_link_libraries(
  sbn
  sbnanalysis_Processor
  #sbnanalysis_Event
  jsoncpp
  ${LARSIM_BASE_LIBRARY}
  ${LARSIM_BASE_DICT}
  ${ROOT_LIBRARIES}
)

install(TARGETS sbnanalysis_Event DESTINATION lib)
install(TARGETS sbnanalysis_Processor DESTINATION lib)
install(TARGETS sbn DESTINATION bin)

if (DEFINED SBNANA_COMPILED_PROCESSOR)
  set(STATICEXE SBNProcessor_${SBNANA_COMPILED_PROCESSOR})
  include_directories( ${CMAKE_INSTALL_PREFIX}/lib)
  link_directories(
    ${CMAKE_INSTALL_PREFIX}/lib
  )
  add_executable(${STATICEXE} sbn_static.cxx)
  target_link_libraries(${STATICEXE} sbnanalysis_Processor jsoncpp sbnanalysis_${SBNANA_COMPILED_PROCESSOR})
  install(TARGETS ${STATICEXE} DESTINATION bin)
endif()

