# Declare libraries to build
# Note that the library names (the first argument) should follow the convention
# "sbnanalysis_ThisAnalysis_LibraryName)". Libraries can be built from multiple
# source files.
add_library(sbnanalysis_SBNOscReco_Histograms SHARED Derived.cc CosmicHisto.cc InteractionHisto.cc TrackHisto.cc DynamicSelector.cc Profile.cc HistoList.cc CRTHisto.cc TrueParticleHisto.cc)

target_link_libraries(sbnanalysis_SBNOscReco_Histograms
  sbnanalysis_SBNOscReco_Reco
  ubana_BasicTool_GeoAlgo
  sbnanalysis_SBNOscReco_classes
  sbnanalysis_SBNOscReco_RecoUtils
  uscript_lib
)

# Install all libraries to the main "lib" directory
install(TARGETS sbnanalysis_SBNOscReco_Histograms DESTINATION lib)


