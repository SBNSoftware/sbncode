# ======================================================================
#  larsoft main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larsoft/ups/setup_for_development <-d|-p>
#  cmake [-DCMAKE_INSTALL_PREFIX=/install/path]
#        -DCMAKE_BUILD_TYPE=$CETPKG_TYPE
#        $CETPKG_SOURCE
#  make
#  make test
#  make install
#  make package (builds distribution tarfile)
# ======================================================================


cmake_minimum_required(VERSION 3.19 FATAL_ERROR)

project(sbncode VERSION 09.34.00 LANGUAGES CXX)

message(STATUS "\n\n  ==========================   ${PROJECT_NAME}   ==========================")

# cetbuildtools contains our cmake modules
find_package(cetbuildtools REQUIRED)

list(APPEND CMAKE_MODULE_PATH $ENV{CANVAS_ROOT_IO_DIR}/Modules $ENV{ART_DIR}/Modules )

include(CetCMakeEnv)
cet_cmake_env()

# The specification of the C++ standard is required because CET macros
# set it to C++98 by default, and then override it depending on the
# compiler version, as determined by the compiler qualifier.
# Removing that qualifier, we need to supply that information manually,
# and maintain it too.
cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  ALLOW_DEPRECATIONS
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs -Wno-undefined-var-template
)

cet_report_compiler_flags()

# these are minimum required versions, not the actual product versions
find_ups_product(art)
find_ups_product(messagefacility)
find_ups_product(art_root_io)
find_ups_product(gsl)
find_ups_product(nusimdata)
find_ups_product(nug4)
find_ups_product(nugen)
find_ups_product(nurandom)
find_ups_product(nutools)
find_ups_product(genie)
find_ups_product(lardataobj)
find_ups_product(larcoreobj)
find_ups_product(larcorealg)
find_ups_product(lardataalg)
find_ups_product(larcore)
find_ups_product(lardata)
find_ups_product(larevt)
find_ups_product(larsim)
find_ups_product(larreco)
find_ups_product(larpandora)
find_ups_product(larsoft v09_00_00)
find_ups_product(artdaq_core v3_06_01)
find_ups_product(ifbeam)
find_ups_product(ifbeam_art)
find_ups_product(libwda) 
find_ups_product(sbnobj v09_10_00)
find_ups_product(sbndata v01_00)
find_ups_product(ifdhc)
find_ups_product(ifdh_art)
find_ups_product(log4cpp)
find_ups_boost()
find_ups_root()
#find_ups_product(sbndcode v06_67_00)
#find_ups_product(uboonecode v06_67_00)
find_ups_product( sbndaq_artdaq_core )

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

add_subdirectory(sbncode)
add_subdirectory(ups)
add_subdirectory(fcl)

# packaging utility
include(UseCPack)
