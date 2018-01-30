# Define the subproject name
project(LeeTruthSelection)

# Directories where includes are found (in addition to those defined in the
# parent CMakeLists.txt, one directory up).
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Declare libraries to build
# Note that the library names (the first argument) should follow the convention
# "sbnanalysis_ThisAnalysis_LibraryName)". Libraries can be built from multiple
# source files.
add_library(sbnanalysis_LeeTruthSelection_TSSelection SHARED TSSelection.cxx)
add_library(sbnanalysis_LeeTruthSelection_TSProcessor SHARED TSProcessor.cxx)
add_library(sbnanalysis_LeeTruthSelection_TSUtil SHARED TSUtil.cxx)

# Define library linking
target_link_libraries(sbnanalysis_LeeTruthSelection_TSSelection sbnanalysis_LeeTruthSelection_TSProcessor sbnanalysis_LeeTruthSelection_TSUtil)

# Install all libraries to the main "lib" directory
install(TARGETS sbnanalysis_LeeTruthSelection_TSSelection DESTINATION lib)
install(TARGETS sbnanalysis_LeeTruthSelection_TSProcessor DESTINATION lib)
install(TARGETS sbnanalysis_LeeTruthSelection_TSUtil DESTINATION lib)

