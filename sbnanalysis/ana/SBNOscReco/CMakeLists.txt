# Define the subproject name
project(SBNOscReco)

# Directories where includes are found (in addition to those defined in the
# parent CMakeLists.txt, one directory up).
include_directories(${CMAKE_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories($ENV{LARCOREALG_INC})
include_directories($ENV{LARSIM_INC})
include_directories($ENV{NUG4_INC})

# build sub-dirs
add_subdirectory(Data)
add_subdirectory(NumuReco)
add_subdirectory(PreSelection)
add_subdirectory(PostProcess)
add_subdirectory(RecoUtils)
add_subdirectory(CosmicIDAlgs)
add_subdirectory(LArReco)
add_subdirectory(OpHitFinder)
add_subdirectory(DataProductViewers)
add_subdirectory(Histograms)
add_subdirectory(TriggerEmulator)
add_subdirectory(fcl)
add_subdirectory(files)
