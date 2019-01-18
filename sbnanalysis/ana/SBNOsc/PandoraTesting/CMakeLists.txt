# Define the subproject name
project(SBNOsc)

# Declare libraries to build
# Note that the library names (the first argument) should follow the convention
# "sbnanalysis_ThisAnalysis_LibraryName)". Libraries can be built from multiple
# source files.
add_library(sbnanalysis_SBNOsc_PandoraViewer SHARED PandoraViewer.cxx)

# Define library linking
target_link_libraries(sbnanalysis_SBNOsc_PandoraViewer
  sbnanalysis_SBNOsc_Utilities
  sbnanalysis_SBNOsc_classes
  ubana_BasicTool_GeoAlgo
)

# Install all libraries to the main "lib" directory
install(TARGETS sbnanalysis_SBNOsc_PandoraViewer DESTINATION lib)

