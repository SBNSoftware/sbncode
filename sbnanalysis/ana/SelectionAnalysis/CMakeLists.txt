# Define the subproject name
project(SelectionAnalysis)

# Directories where includes are found (in addition to those defined in the
# parent CMakeLists.txt, one directory up).
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# FOR LARSOFT CMAKE NEEDS
include_directories($ENV{LARCOREALG_INC})
include_directories($ENV{LARSIM_INC})

# FOR LARSOFT CMAKE NEEDS
include_directories($ENV{LARCOREALG_INC})
include_directories($ENV{LARSIM_INC})

# Declare libraries to build
# Note that the library names (the first argument) should follow the convention
# "sbnanalysis_ThisAnalysis_LibraryName)". Libraries can be built from multiple
add_library(sbnanalysis_SelectionAnalysis_AnalysisMain SHARED AnalysisMain.cxx)

# Add in ROOT dictionary
ROOT_GENERATE_DICTIONARY(SelectionAnalysis_dict AnalysisMain.h LINKDEF linkdef.h)
add_library(sbnanalysis_SelectionAnalysis_classes SHARED SelectionAnalysis_dict.cxx)
add_custom_command(TARGET sbnanalysis_SelectionAnalysis_classes POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
 libSelectionAnalysis_dict_rdict.pcm ${CMAKE_INSTALL_PREFIX}/lib)
target_link_libraries(sbnanalysis_SelectionAnalysis_AnalysisMain
  sbnanalysis_SelectionAnalysis_classes
  ubana_BasicTool_GeoAlgo  
)
##add_custom_command(TARGET sbnanalysis_SelectionAnalysis_classes POST_BUILD
##  COMMAND ${CMAKE_COMMAND} -E copy_if_different
## libSelectionAnalysis_dict_rdict.pcm ${CMAKE_INSTALL_PREFIX}/lib)

target_link_libraries(sbnanalysis_SelectionAnalysis_AnalysisMain
  sbnanalysis_SelectionAnalysis_classes
  ubana_BasicTool_GeoAlgo  
)

# Install all libraries to the main "lib" directory
install(TARGETS sbnanalysis_SelectionAnalysis_AnalysisMain DESTINATION lib)
install(TARGETS sbnanalysis_SelectionAnalysis_classes DESTINATION lib)

add_subdirectory(config)
add_library(sbnanalysis_SelectionAnalysis_TestDuringMerge SHARED TestDuringMerge.cxx)

# Define library linking
target_link_libraries(sbnanalysis_SelectionAnalysis_TestDuringMerge sbnanalysis_SelectionTool_ObjectBuilder)

# Install all libraries to the main "lib" directory
install(TARGETS sbnanalysis_SelectionAnalysis_TestDuringMerge DESTINATION lib)
