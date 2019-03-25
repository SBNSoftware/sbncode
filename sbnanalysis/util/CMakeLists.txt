include_directories( ${ROOT_INCLUDE_DIRS} )
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
include_directories($ENV{UBOONECODE_INC})
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
  $ENV{LARSIM_LIB}
  $ENV{SBNDCODE_LIB}
  $ENV{UBCORE_LIB}
  $ENV{ICARUSCODE_LIB}
  $ENV{NUTOOLS_LIB}
  $ENV{UBOONECODE_LIB}
  ${CMAKE_INSTALL_DIR}/lib
)

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_library(sbnanalysis_util_Interaction Interaction.cxx)
target_link_libraries(sbnanalysis_util_Interaction ${ROOT_LIBRARIES})

add_library(sbnanalysis_util_RecoUtil RecoUtil.cxx)
target_link_libraries(
  sbnanalysis_util_RecoUtil 
  #sbnanalysis_core_ProviderManager
  #larcorealg_Geometry
  #lardataalg_DetectorInfo
  #larsim_MCCheater_ParticleInventory
  #larsim_MCCheater_BackTracker
  #icaruscode_Geometry
  #ubcore_Geometry
  #sbnd_Geometry
  #boost_system
  #boost_filesystem
  ${ROOT_LIBRARIES}
  )

install(TARGETS sbnanalysis_util_Interaction DESTINATION lib)
install(TARGETS sbnanalysis_util_RecoUtil DESTINATION lib)

