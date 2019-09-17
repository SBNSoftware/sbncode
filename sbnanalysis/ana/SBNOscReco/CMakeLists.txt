# Define the subproject name
project(SBNOscReco)


# build sub-dirs
add_subdirectory(RecoUtils)
add_subdirectory(CosmicIDAlgs)
add_subdirectory(LArReco)
add_subdirectory(OpHitFinder)
add_subdirectory(fcl)

# Directories where includes are found (in addition to those defined in the
# parent CMakeLists.txt, one directory up).
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories($ENV{LARCOREALG_INC})
include_directories($ENV{LARSIM_INC})
include_directories($ENV{NUG4_INC})

# Declare libraries to build
# Note that the library names (the first argument) should follow the convention
# "sbnanalysis_ThisAnalysis_LibraryName)". Libraries can be built from multiple
# source files.
add_library(sbnanalysis_SBNOscReco_NumuRecoSelection SHARED NumuRecoSelection.cxx)
add_library(sbnanalysis_SBNOscReco_GetFV SHARED GetFV.cc)
add_library(sbnanalysis_SBNOscReco_OpDetWaveformMaker SHARED OpDetWaveformMaker.cc)
add_library(sbnanalysis_SBNOscReco_OpSimHitPrinter SHARED OpSimHitPrinter.cc)
add_library(sbnanalysis_SBNOscReco_SimPhotonPrinter SHARED SimPhotonPrinter.cc)
add_library(sbnanalysis_SBNOscReco_OpDetBacktrackerPrinter SHARED OpDetBacktrackerPrinter.cc)
add_library(sbnanalysis_SBNOscReco_MCParticleTreePrinter SHARED MCParticleTreePrinter.cc)
add_library(sbnanalysis_SBNOscReco_PandoraIDPrinter SHARED PandoraIDPrinter.cc)
add_library(sbnanalysis_SBNOscReco_CRTSimHitViewer SHARED CRTSimHitViewer.cc)

# Add in ROOT dictionary
ROOT_GENERATE_DICTIONARY(SBNOscReco_dict NumuRecoSelection.h LINKDEF linkdef.h)
add_library(sbnanalysis_SBNOscReco_classes SHARED SBNOscReco_dict.cxx)
add_custom_command(TARGET sbnanalysis_SBNOscReco_classes POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
  libSBNOscReco_dict_rdict.pcm ${CMAKE_INSTALL_PREFIX}/lib)

target_link_libraries(sbnanalysis_SBNOscReco_NumuRecoSelection
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

target_link_libraries(sbnanalysis_SBNOscReco_GetFV
  sbnanalysis_SBNOsc_Utilities
  sbnanalysis_SBNOscReco_classes
  ubana_BasicTool_GeoAlgo
)

target_link_libraries(sbnanalysis_SBNOscReco_OpDetWaveformMaker
  sbnanalysis_SBNOsc_Utilities
  sbnanalysis_SBNOscReco_classes
  ubana_BasicTool_GeoAlgo
)

target_link_libraries(sbnanalysis_SBNOscReco_OpSimHitPrinter
  sbnanalysis_SBNOsc_Utilities
  sbnanalysis_SBNOscReco_classes
  ubana_BasicTool_GeoAlgo
  sbnanalysis_SBNOscReco_opHitFinderSBND
)

target_link_libraries(sbnanalysis_SBNOscReco_SimPhotonPrinter
  sbnanalysis_SBNOsc_Utilities
  sbnanalysis_SBNOscReco_classes
  ubana_BasicTool_GeoAlgo
)

target_link_libraries(sbnanalysis_SBNOscReco_OpDetBacktrackerPrinter
  sbnanalysis_SBNOsc_Utilities
  sbnanalysis_SBNOscReco_classes
  ubana_BasicTool_GeoAlgo
)

target_link_libraries(sbnanalysis_SBNOscReco_MCParticleTreePrinter
  sbnanalysis_SBNOsc_Utilities
  sbnanalysis_SBNOscReco_classes
  ubana_BasicTool_GeoAlgo
)

target_link_libraries(sbnanalysis_SBNOscReco_PandoraIDPrinter
  sbnanalysis_SBNOsc_Utilities
  sbnanalysis_SBNOscReco_classes
  ubana_BasicTool_GeoAlgo
)

target_link_libraries(sbnanalysis_SBNOscReco_CRTSimHitViewer
  sbnanalysis_SBNOsc_Utilities
  sbnanalysis_SBNOscReco_classes
  ubana_BasicTool_GeoAlgo
  icaruscode_CRTProducts
)

# Install all libraries to the main "lib" directory
install(TARGETS sbnanalysis_SBNOscReco_NumuRecoSelection DESTINATION lib)
install(TARGETS sbnanalysis_SBNOscReco_GetFV DESTINATION lib)
install(TARGETS sbnanalysis_SBNOscReco_OpDetWaveformMaker DESTINATION lib)
install(TARGETS sbnanalysis_SBNOscReco_OpSimHitPrinter DESTINATION lib)
install(TARGETS sbnanalysis_SBNOscReco_SimPhotonPrinter DESTINATION lib)
install(TARGETS sbnanalysis_SBNOscReco_OpDetBacktrackerPrinter DESTINATION lib)
install(TARGETS sbnanalysis_SBNOscReco_MCParticleTreePrinter DESTINATION lib)
install(TARGETS sbnanalysis_SBNOscReco_PandoraIDPrinter DESTINATION lib)
install(TARGETS sbnanalysis_SBNOscReco_CRTSimHitViewer DESTINATION lib)
install(TARGETS sbnanalysis_SBNOscReco_classes DESTINATION lib)

