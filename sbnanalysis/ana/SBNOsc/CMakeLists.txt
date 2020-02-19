# Define the subproject name
project(SBNOsc)

# Directories where includes are found (in addition to those defined in the
# parent CMakeLists.txt, one directory up).
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories($ENV{LARCOREALG_INC})
include_directories($ENV{LARSIM_INC})

# Declare libraries to build
# Note that the library names (the first argument) should follow the convention
# "sbnanalysis_ThisAnalysis_LibraryName)". Libraries can be built from multiple
# source files.
add_library(sbnanalysis_SBNOsc_Utilities SHARED Utilities.cxx)
add_library(sbnanalysis_SBNOsc_NueSelection SHARED NueSelection.cxx)
add_library(sbnanalysis_SBNOsc_NumuSelection SHARED NumuSelection.cxx)
add_library(sbnanalysis_SBNOsc_Covariance SHARED Covariance.cxx)
add_library(sbnanalysis_SBNOsc_Chi2Sensitivity SHARED Chi2Sensitivity.cxx)

# Add in ROOT dictionary
ROOT_GENERATE_DICTIONARY(SBNOsc_dict NumuSelection.h LINKDEF linkdef.h)
add_library(sbnanalysis_SBNOsc_classes SHARED SBNOsc_dict.cxx)
add_custom_command(TARGET sbnanalysis_SBNOsc_classes POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
  libSBNOsc_dict_rdict.pcm ${CMAKE_INSTALL_PREFIX}/lib)

# Define library linking
target_link_libraries(sbnanalysis_SBNOsc_NueSelection sbnanalysis_SBNOsc_Utilities)
target_link_libraries(sbnanalysis_SBNOsc_NumuSelection
  sbnanalysis_SBNOsc_Utilities
  sbnanalysis_SBNOsc_classes
  ubana_BasicTool_GeoAlgo
)
target_link_libraries(sbnanalysis_SBNOsc_Chi2Sensitivity sbnanalysis_SBNOsc_Covariance)

# Install all libraries to the main "lib" directory
install(TARGETS sbnanalysis_SBNOsc_Utilities DESTINATION lib)
install(TARGETS sbnanalysis_SBNOsc_NueSelection DESTINATION lib)
install(TARGETS sbnanalysis_SBNOsc_NumuSelection DESTINATION lib)
install(TARGETS sbnanalysis_SBNOsc_Covariance DESTINATION lib)
install(TARGETS sbnanalysis_SBNOsc_Chi2Sensitivity DESTINATION lib)

install(TARGETS sbnanalysis_SBNOsc_classes DESTINATION lib)

add_subdirectory(config)


# build the pandora viewer
add_library(sbnanalysis_SBNOsc_PandoraViewer SHARED PandoraViewer.cxx)

target_link_libraries(sbnanalysis_SBNOsc_PandoraViewer
  sbnanalysis_SBNOsc_Utilities
  ubana_BasicTool_GeoAlgo
  sbndcode_CRT
  sbndcode_CRTProducts
  sbndcode_CRTUtils
  sbndcode_CRT_CRTProducts_dict
)

install(TARGETS sbnanalysis_SBNOsc_PandoraViewer DESTINATION lib)
