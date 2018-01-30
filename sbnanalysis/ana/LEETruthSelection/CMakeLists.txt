# Define the subproject name
project(LEETruthSelection)

# Directories where includes are found (in addition to those defined in the
# parent CMakeLists.txt, one directory up).
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Declare libraries to build
# Note that the library names (the first argument) should follow the convention
# "sbnanalysis_ThisAnalysis_LibraryName)". Libraries can be built from multiple
# source files.
add_library(sbnanalysis_LEETruthSelection_LEESelection SHARED LEESelection.cxx)
add_library(sbnanalysis_LEETruthSelection_Config SHARED Config.cxx)
add_library(sbnanalysis_LEETruthSelection_Cuts SHARED Cuts.cxx)
add_library(sbnanalysis_LEETruthSelection_Util SHARED Util.cxx)

# Define library linking
target_link_libraries(
  sbnanalysis_LEETruthSelection_LEESelection
  sbnanalysis_Processor
  sbnanalysis_Event
  sbnanalysis_LEETruthSelection_Config
  sbnanalysis_LEETruthSelection_Cuts
  sbnanalysis_LEETruthSelection_Util
)

# Install all libraries to the main "lib" directory
install(TARGETS sbnanalysis_LEETruthSelection_LEESelection DESTINATION lib)
install(TARGETS sbnanalysis_LEETruthSelection_Config DESTINATION lib)
install(TARGETS sbnanalysis_LEETruthSelection_Cuts DESTINATION lib)
install(TARGETS sbnanalysis_LEETruthSelection_Util DESTINATION lib)

