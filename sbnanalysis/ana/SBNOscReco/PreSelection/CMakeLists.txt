# Declare libraries to build
# Note that the library names (the first argument) should follow the convention
# "sbnanalysis_ThisAnalysis_LibraryName)". Libraries can be built from multiple
# source files.
add_library(sbnanalysis_SBNOscReco_NumuReco SHARED NumuReco.cxx)

target_link_libraries(sbnanalysis_SBNOscReco_NumuReco
  sbnanalysis_SBNOscReco_Reco
  sbnanalysis_SBNOscReco_Histograms
  sbnanalysis_SBNOsc_Utilities
  sbnanalysis_SBNOscReco_classes
  sbnanalysis_SBNOscReco_RecoUtils
  sbnanalysis_SBNOscReco_LArReco
  ubana_BasicTool_GeoAlgo
  sbncode_CosmicIdAlgs
  sbnanalysis_SBNOscReco_opHitFinderSBND
  icaruscode_CRT	
  icaruscode_CRTUtils
  icaruscode_CRTData
  icaruscode_CRTProducts
  icaruscode_CRT_CRTProducts_dict
)

# Install all libraries to the main "lib" directory
install(TARGETS sbnanalysis_SBNOscReco_NumuReco DESTINATION lib)


