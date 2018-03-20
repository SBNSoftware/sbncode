# Define the subproject name
project(SBNProposalAnalysis)

# Directories where includes are found (in addition to those defined in the
# parent CMakeLists.txt, one directory up).
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Declare libraries to build
# Note that the library names (the first argument) should follow the convention
# "sbnanalysis_ThisAnalysis_LibraryName)". Libraries can be built from multiple
# source files.
add_library(sbnanalysis_SBNProposalAnalysis_SBNProposalTools SHARED SBNProposalTools.cxx)
add_library(sbnanalysis_SBNProposalAnalysis_SBNProposalSelection SHARED SBNProposalSelection.cxx)

# Define library linking
target_link_libraries(sbnanalysis_SBNProposalAnalysis_SBNProposalSelection sbnanalysis_SBNProposalAnalysis_SBNProposalTools)

# Install all libraries to the main "lib" directory
install(TARGETS sbnanalysis_SBNProposalAnalysis_SBNProposalTools DESTINATION lib)
install(TARGETS sbnanalysis_SBNProposalAnalysis_SBNProposalSelection DESTINATION lib)

