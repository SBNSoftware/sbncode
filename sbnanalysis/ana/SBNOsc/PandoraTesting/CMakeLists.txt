# Define the subproject name
project(SBNOsc)

# Declare libraries to build
# Note that the library names (the first argument) should follow the convention
# "sbnanalysis_ThisAnalysis_LibraryName)". Libraries can be built from multiple
# source files.
add_library(sbnanalysis_SBNOsc_PandoraViewer SHARED PandoraViewer.cxx)
add_library(sbnanalysis_SBNOsc_PandoraTesting SHARED PandoraTesting.cxx)

# Define library linking
target_link_libraries(sbnanalysis_SBNOsc_PandoraViewer
  sbnanalysis_SBNOsc_Utilities
  sbnanalysis_PandoraTesting_classes
  ubana_BasicTool_GeoAlgo
  sbndcode_CRT
  sbndcode_CRTProducts
  sbndcode_CRTUtils
  sbndcode_CRT_CRTProducts_dict
)

target_link_libraries(sbnanalysis_SBNOsc_PandoraTesting
  sbnanalysis_SBNOsc_Utilities
  sbnanalysis_PandoraTesting_classes
  ubana_BasicTool_GeoAlgo
  sbndcode_CRT
  sbndcode_CRTProducts
  sbndcode_CRTUtils
  sbndcode_CRT_CRTProducts_dict
)

# Add in ROOT dictionary
ROOT_GENERATE_DICTIONARY(PandoraTesting_dict PandoraTesting.h LINKDEF linkdef.h)
add_library(sbnanalysis_PandoraTesting_classes SHARED PandoraTesting_dict.cxx)
add_custom_command(TARGET sbnanalysis_PandoraTesting_classes POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
  libPandoraTesting_dict_rdict.pcm ${CMAKE_INSTALL_PREFIX}/lib)

# Install all libraries to the main "lib" directory
install(TARGETS sbnanalysis_SBNOsc_PandoraViewer DESTINATION lib)
install(TARGETS sbnanalysis_SBNOsc_PandoraTesting DESTINATION lib)
install(TARGETS sbnanalysis_PandoraTesting_classes DESTINATION lib)

