# Define the subproject name
project(AutoVetoAnalysis)

# Directories where includes are found (in addition to those defined in the
# parent CMakeLists.txt, one directory up).
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
#include_directories( $ENV{ICARUSCODE_INC} )

#link_directories(
#    $ENV{ICARUSCODE_LIB}
#)

# Declare libraries to build
# Note that the library names (the first argument) should follow the convention
# "sbnanalysis_ThisAnalysis_LibraryName)". Libraries can be built from multiple
# source files.
add_library(sbnanalysis_AutoVetoAnalysis_AutoVetoTools SHARED AutoVetoTools.cxx)
add_library(sbnanalysis_AutoVetoAnalysis_AutoVetoSelection SHARED AutoVetoSelection.cxx)

# Define library linking
target_link_libraries(
  sbnanalysis_AutoVetoAnalysis_AutoVetoSelection
  sbnanalysis_AutoVetoAnalysis_AutoVetoTools
  sbnanalysis_ProviderManager
  #icaruscode_CRT
  #icaruscode_CRTData
  #icaruscode_CRTProducts
)

# Install all libraries to the main "lib" directory
install(TARGETS sbnanalysis_AutoVetoAnalysis_AutoVetoTools DESTINATION lib)
install(TARGETS sbnanalysis_AutoVetoAnalysis_AutoVetoSelection DESTINATION lib)

