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


cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

find_package(cetmodules 3.20.00 REQUIRED)
project(sbncode VERSION 09.93.01 LANGUAGES CXX)

message(STATUS "\n\n  ==========================   ${PROJECT_NAME}   ==========================")


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
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
)
set(fcl_dir caf)
cet_report_compiler_flags(REPORT_THRESHOLD VERBOSE)

# these are minimum required versions, not the actual product versions
find_package( art REQUIRED )
find_package( messagefacility REQUIRED )
find_package( art_root_io REQUIRED )
find_package( nusimdata REQUIRED )
find_package( nug4 REQUIRED )
find_package( nugen REQUIRED )
find_package( nurandom REQUIRED )
find_package( nutools REQUIRED )
find_package( systematicstools REQUIRED )
find_package( nusystematics REQUIRED )
find_package( GENIE REQUIRED )
find_package( lardataobj REQUIRED )
find_package( larcoreobj REQUIRED )
find_package( larcorealg REQUIRED )
find_package( lardataalg REQUIRED )
find_package( larcore REQUIRED )
find_package( lardata REQUIRED )
find_package( larevt REQUIRED )
find_package( larsim REQUIRED )
find_package( larreco REQUIRED )
find_package( larrecodnn REQUIRED )
find_package( larpandora REQUIRED )
find_package( larpandoracontent REQUIRED )
find_package( artdaq_core EXPORT REQUIRED )
find_package( sbnobj REQUIRED )
find_package( sbnanaobj REQUIRED )
find_package( sbndaq_artdaq_core EXPORT REQUIRED )
find_package( ifdh_art REQUIRED )
find_package( ifbeam REQUIRED )
find_package( libwda REQUIRED )
find_package( dk2nudata REQUIRED )
find_package( dk2nugenie EXPORT REQUIRED )
find_package( log4cpp REQUIRED )
find_package( LibXml2 REQUIRED )
find_package( CLHEP REQUIRED )
find_package( ROOT REQUIRED )
find_package( Geant4 REQUIRED )
find_package( Boost COMPONENTS system filesystem REQUIRED )

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

add_subdirectory(sbncode)
add_subdirectory(fcl)
add_subdirectory(fclmodule)


#generate cmake config files
cet_cmake_config()
