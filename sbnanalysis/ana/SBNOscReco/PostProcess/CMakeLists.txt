# Declare libraries to build
# Note that the library names (the first argument) should follow the convention
# "sbnanalysis_ThisAnalysis_LibraryName)". Libraries can be built from multiple
# source files.
add_library(sbnanalysis_SBNOscReco_Selection SHARED Selection.cc Cuts.cc Histograms.cc)

target_link_libraries(sbnanalysis_SBNOscReco_Selection
  sbnanalysis_SBNOsc_Utilities
  ubana_BasicTool_GeoAlgo
  sbnanalysis_SBNOscReco_classes
)

# Install all libraries to the main "lib" directory
install(TARGETS sbnanalysis_SBNOscReco_Selection DESTINATION lib)


