include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

include_directories(${ROOT_INCLUDE_DIRS})
include_directories($ENV{ART_INC})
include_directories($ENV{BOOST_INC})
include_directories($ENV{NUTOOLS_INC})
include_directories($ENV{FHICLCPP_INC})
include_directories($ENV{MESSAGEFACILITY_INC})
include_directories($ENV{CANVAS_INC})
include_directories($ENV{CANVAS_ROOT_IO_INC})
include_directories($ENV{CETLIB_INC})
include_directories($ENV{CETLIB_EXCEPT_INC})
include_directories($ENV{GALLERY_INC})
include_directories($ENV{CLHEP_INC})
include_directories($ENV{LARCORE_INC})
include_directories($ENV{LARCOREOBJ_INC})
include_directories($ENV{LARRECO_INC})
include_directories($ENV{LARCORE_INC})
include_directories($ENV{LARCOREALG_INC})
include_directories($ENV{LARDATA_INC})
include_directories($ENV{LARDATAOBJ_INC})
include_directories($ENV{LARDATAALG_INC})
include_directories($ENV{NUSIMDATA_INC})
include_directories($ENV{LARSIM_INC})
include_directories($ENV{SBNDCODE_INC})
include_directories($ENV{UBCORE_INC})
include_directories($ENV{ICARUSCODE_INC})
include_directories($ENV{MESSAGEFACILITY_INC})

link_directories(
  $ENV{ART_LIB}
  $ENV{BOOST_LIB}
  $ENV{CANVAS_LIB}
  $ENV{CETLIB_LIB}
  $ENV{GALLERY_LIB}
  $ENV{NUSIMDATA_LIB}
  $ENV{LARCORE_LIB}
  $ENV{LARCOREOBJ_LIB}
  $ENV{LARCOREALG_LIB}
  $ENV{LARRECO_LIB}
  $ENV{LARDATA_LIB}
  $ENV{LARDATAOBJ_LIB}
  $ENV{LARDATAALG_LIB}
  $ENV{FHICLCPP_LIB}
  $ENV{MESSAGEFACILITY_LIB}
  $ENV{LARDATAALG_LIB}
  $ENV{JSONCPP_LIB}
  $ENV{MESSAGEFACILITY_LIB}
  $ENV{LARSIM_LIB}
  $ENV{SBNDCODE_LIB}
  $ENV{UBCORE_LIB}
  $ENV{ICARUSCODE_LIB}
  $ENV{NUTOOLS_LIB}
)

ROOT_GENERATE_DICTIONARY(sbn_eventdict Event.hh SubRun.hh LINKDEF linkdef.h)

add_library(sbnanalysis_Event SHARED sbn_eventdict.cxx)
target_link_libraries(
  sbnanalysis_Event
  larsim_EventWeight_Base
  larsim_EventWeight_Base_dict
  larsim_Simulation 
  lardataobj_Simulation
  ${ROOT_LIBRARIES}
)

add_library(
  sbnanalysis_Processor SHARED
  ProcessorBase.cxx
  SelectionBase.cxx
  ProcessorBlock.cxx
  Loader.cxx
  PostProcessorBase.cxx
)

add_library(
  sbnanalysis_ProviderManager SHARED
  ProviderManager.cxx
)

target_link_libraries(
  sbnanalysis_ProviderManager
  larcorealg_Geometry
  lardataalg_DetectorInfo
  larsim_MCCheater_ParticleInventory
  larsim_MCCheater_BackTracker
  sbnanalysis_ServiceManager SHARED
  ServiceManager.cxx
  sbnanalysis_ProviderManager SHARED
  ProviderManager.cxx
)

target_link_libraries(
  sbnanalysis_ProviderManager
  larcorealg_Geometry
  lardataalg_DetectorInfo
  larsim_MCCheater_ParticleInventory
  larsim_MCCheater_BackTracker
  icaruscode_Geometry
  ubcore_Geometry
  sbnd_Geometry
  boost_system
  boost_filesystem
)

target_link_libraries(
  sbnanalysis_Processor
  sbnanalysis_Event
  sbnanalysis_util_Interaction
  sbnanalysis_util_RecoUtil
  sbnanalysis_ProviderManager
  sbnanalysis_ServiceManager
  canvas
  cetlib
  gallery
  nusimdata_SimulationBase
  larcoreobj_SummaryData
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  lardataobj_RecoBase
  lardataobj_AnalysisBase
  lardataobj_MCBase
  larsim_EventWeight_Base
  larsim_EventWeight_Base_dict
  larsim_Simulation 
  lardataobj_Simulation
  lardataalg_DetectorInfo
  larreco_Calorimetry
  larreco_RecoAlg
  larreco_RecoAlg_ClusterRecoUtil
  larreco_ClusterFinder
  sbndcode_RecoUtils
  ${CLHEP_Geometry}
  art_Persistency_Common canvas
  art_Utilities canvas
  ${ART_FRAMEWORK_CORE}                                  
  ${ART_FRAMEWORK_PRINCIPAL}                             
  ${ART_FRAMEWORK_SERVICES_REGISTRY}                     
  ${ART_FRAMEWORK_SERVICES_BASIC}                        
  ${ART_FRAMEWORK_SERVICES_OPTIONAL}                     
  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
  ${MF_MESSAGELOGGER}
  ${ROOT_LIBRARIES}
)

file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/lib)
add_custom_command(TARGET sbnanalysis_Event POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
  libsbn_eventdict_rdict.pcm ${CMAKE_INSTALL_PREFIX}/lib)

file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/fcl)
add_custom_command(TARGET sbnanalysis_ProviderManager POST_BUILD
add_custom_command(TARGET sbnanalysis_ServiceManager POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_directory
  ${CMAKE_CURRENT_SOURCE_DIR}/fcl ${CMAKE_INSTALL_PREFIX}/fcl)

add_executable(sbn sbn.cxx)
target_link_libraries(
  sbn
  sbnanalysis_Processor
  ${LARSIM_BASE_LIBRARY}
  ${LARSIM_BASE_DICT}
  fhiclcpp
  nusimdata_SimulationBase
  larcoreobj_SummaryData
  larcorealg_Geometry
  larcore_Geometry_Geometry_service
  lardataobj_RecoBase
  lardataobj_MCBase
  larsim_EventWeight_Base
  larsim_EventWeight_Base_dict
  larsim_Simulation 
  lardataobj_Simulation
  lardataalg_DetectorInfo
  larreco_Calorimetry
  larreco_RecoAlg
  larreco_RecoAlg_ClusterRecoUtil
  larreco_ClusterFinder
  sbndcode_RecoUtils
  art_Persistency_Common canvas
  art_Utilities canvas
  ${CLHEP_Geometry}
  ${ART_FRAMEWORK_CORE}                                  
  ${ART_FRAMEWORK_PRINCIPAL}                             
  ${ART_FRAMEWORK_SERVICES_REGISTRY}                     
  ${ART_FRAMEWORK_SERVICES_BASIC}                        
  ${ART_FRAMEWORK_SERVICES_OPTIONAL}                     
  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
  ${MF_MESSAGELOGGER}
  jsoncpp
  ${LARSIM_BASE_LIBRARY}
  ${LARSIM_BASE_DICT}
  larsim_Simulation 
  lardataobj_Simulation
  lardataalg_DetectorInfo
  larreco_Calorimetry
  larreco_RecoAlg
  larreco_RecoAlg_ClusterRecoUtil
  larreco_ClusterFinder
  sbndcode_RecoUtils
  art_Persistency_Common canvas
  art_Utilities canvas
  ${CLHEP_Geometry}
  ${ART_FRAMEWORK_CORE}                                  
  ${ART_FRAMEWORK_PRINCIPAL}                             
  ${ART_FRAMEWORK_SERVICES_REGISTRY}                     
  ${ART_FRAMEWORK_SERVICES_BASIC}                        
  ${ART_FRAMEWORK_SERVICES_OPTIONAL}                     
  ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
  ${MF_MESSAGELOGGER}
  ${ROOT_LIBRARIES}
)

add_executable(sbn-postprocess sbn-postprocess.cxx)
target_link_libraries(
  sbn-postprocess
  sbnanalysis_Processor
  ${LARSIM_BASE_LIBRARY}
  ${LARSIM_BASE_DICT}
  fhiclcpp
  larsim_EventWeight_Base
  larsim_EventWeight_Base_dict
  larsim_Simulation 
  jsoncpp
  ${LARSIM_BASE_LIBRARY}
  ${LARSIM_BASE_DICT}
  larsim_Simulation 
  ${ROOT_LIBRARIES}
)

install(TARGETS sbnanalysis_Event DESTINATION lib)
install(TARGETS sbnanalysis_Processor DESTINATION lib)
install(TARGETS sbnanalysis_ProviderManager DESTINATION lib)
install(TARGETS sbnanalysis_ServiceManager DESTINATION lib)
install(TARGETS sbn DESTINATION bin)
install(TARGETS sbn-postprocess DESTINATION bin)

if (DEFINED SBNANA_COMPILED_PROCESSOR)
  set(STATICEXE SBNProcessor_${SBNANA_COMPILED_PROCESSOR})
  include_directories( ${CMAKE_INSTALL_PREFIX}/lib)
  link_directories(
    ${CMAKE_INSTALL_PREFIX}/lib
  )
  add_executable(${STATICEXE} sbn_static.cxx)
  target_link_libraries(${STATICEXE} sbnanalysis_Processor fhiclcpp sbnanalysis_${SBNANA_COMPILED_PROCESSOR})
  install(TARGETS ${STATICEXE} DESTINATION bin)
endif()
