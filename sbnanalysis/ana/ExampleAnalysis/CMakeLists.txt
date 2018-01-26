# Define the subproject name
project(ExampleAnalysis)

# Directories where includes are found (in addition to those defined in the
# parent CMakeLists.txt, one directory up).
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Declare libraries to build
# Note that the library names (the first argument) should follow the convention
# "sbnanalysis_ThisAnalysis_LibraryName)". Libraries can be built from multiple
# source files.
add_library(sbnanalysis_ExampleAnalysis_ExampleTools SHARED ExampleTools.cxx)
add_library(sbnanalysis_ExampleAnalysis_ExampleSelection SHARED ExampleSelection.cxx)

# Define library linking
target_link_libraries(sbnanalysis_ExampleAnalysis_ExampleSelection sbnanalysis_ExampleAnalysis_ExampleTools)

# Install all libraries to the main "lib" directory
install(TARGETS sbnanalysis_ExampleAnalysis_ExampleTools DESTINATION lib)
install(TARGETS sbnanalysis_ExampleAnalysis_ExampleSelection DESTINATION lib)

