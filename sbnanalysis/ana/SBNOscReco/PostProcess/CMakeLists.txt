# Declare libraries to build
# Note that the library names (the first argument) should follow the convention
# "sbnanalysis_ThisAnalysis_LibraryName)". Libraries can be built from multiple
# source files.
add_library(sbnanalysis_SBNOscReco_Selection SHARED Selection.cc Normalize.cc ROC.cc Histograms.cc Cuts.cc SetEvent.cc)

add_library(sbnanalysis_SBNOscReco_PIDDataMaker SHARED PIDDataMaker.cc Cuts.cc SetEvent.cc)

add_library(sbnanalysis_SBNOscReco_Flatten SHARED Flatten.cc Cuts.cc SetEvent.cc)

target_link_libraries(sbnanalysis_SBNOscReco_Selection
  sbnanalysis_SBNOscReco_TriggerEmulator
  sbnanalysis_SBNOscReco_Histograms
  ubana_BasicTool_GeoAlgo
  sbndcode_GeoWrappers
  sbnanalysis_SBNOscReco_classes
  sbnanalysis_SBNOscReco_RecoUtils
  sbnanalysis_Event
  xgboost
)

target_link_libraries(sbnanalysis_SBNOscReco_Flatten
  sbnanalysis_SBNOscReco_classes
  sbnanalysis_SBNOscReco_TriggerEmulator
  sbnanalysis_SBNOscReco_Histograms
  sbnanalysis_SBNOscReco_RecoUtils
)

target_link_libraries(sbnanalysis_SBNOscReco_PIDDataMaker
  sbnanalysis_SBNOscReco_classes
  sbnanalysis_SBNOscReco_TriggerEmulator
  sbnanalysis_SBNOscReco_Histograms
  sbnanalysis_SBNOscReco_RecoUtils
)

# Install all libraries to the main "lib" directory
install(TARGETS sbnanalysis_SBNOscReco_Selection DESTINATION lib)

install(TARGETS sbnanalysis_SBNOscReco_PIDDataMaker DESTINATION lib)
install(TARGETS sbnanalysis_SBNOscReco_Flatten DESTINATION lib)
