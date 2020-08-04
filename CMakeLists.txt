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


# use cmake 2.8 or later
cmake_minimum_required (VERSION 2.8)

project(sbncode)

message(STATUS "\n\n  ==========================   ${PROJECT_NAME}   ==========================")

# cetbuildtools contains our cmake modules
SET ( CETBUILDTOOLS_VERSION $ENV{CETBUILDTOOLS_VERSION} )
IF (NOT CETBUILDTOOLS_VERSION)
    MESSAGE (FATAL_ERROR "ERROR: setup cetbuildtools to get the cmake modules")
ENDIF()

set(CMAKE_MODULE_PATH $ENV{CANVAS_ROOT_IO_DIR}/Modules
                      $ENV{CETBUILDTOOLS_DIR}/Modules
                      ${CMAKE_MODULE_PATH})

include(CetCMakeEnv)
cet_cmake_env()

set_install_root()

# The specification of the C++ standard is required because CET macros
# set it to C++98 by default, and then override it depending on the
# compiler version, as determined by the compiler qualifier.
# Removing that qualifier, we need to supply that information manually,
# and maintain it too.
cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  ALLOW_DEPRECATIONS
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)

cet_report_compiler_flags()

# these are minimum required versions, not the actual product versions
find_ups_product( larpandora v1_00_00 )
find_ups_product( larsim v1_00_00 )
find_ups_product( lardata v1_00_00 )
find_ups_product( lardataobj )
find_ups_product( lardataalg )
find_ups_product( larevt )
find_ups_product( larcore v1_00_00 )
find_ups_product( larcoreobj v1_00_00 )
find_ups_product( larreco )
find_ups_product( nutools v1_00_00 )
find_ups_product( nug4 )
find_ups_product( nugen )
find_ups_product( nuevdb )
find_ups_product( nurandom v1_00_00 )
find_ups_product( art v1_09_00 )
find_ups_product( art_root_io )
find_ups_product( postgresql v9_1_5 )
find_ups_product( hep_concurrency )
find_ups_product( pandora )
find_ups_product( eigen )

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net)
include(${ROOT_USE_FILE})

add_subdirectory(sbncode)
add_subdirectory(ups)
add_subdirectory(fcl)

# packaging utility
include(UseCPack)

