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
project(sbncode VERSION 10.31.00.01 LANGUAGES CXX)

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
find_package( art REQUIRED EXPORT )
find_package( Eigen3 REQUIRED EXPORT)
find_package( GSL REQUIRED EXPORT )
find_package( messagefacility REQUIRED EXPORT )
find_package( art_root_io REQUIRED EXPORT )
find_package( nusimdata REQUIRED EXPORT )
find_package( nug4 REQUIRED EXPORT )
find_package( nugen REQUIRED EXPORT )
find_package( nurandom REQUIRED EXPORT )
find_package( nutools REQUIRED EXPORT )
find_package( systematicstools REQUIRED EXPORT )
find_package( nusystematics REQUIRED EXPORT )
find_package( GENIE REQUIRED EXPORT )
find_package( lardataobj REQUIRED EXPORT )
find_package( larcoreobj REQUIRED EXPORT )
find_package( larcorealg REQUIRED EXPORT )
find_package( lardataalg REQUIRED EXPORT )
find_package( larcore REQUIRED EXPORT )
find_package( lardata REQUIRED EXPORT )
find_package( larevt REQUIRED EXPORT )
find_package( larsim REQUIRED EXPORT )
find_package( larreco REQUIRED EXPORT )
find_package( larrecodnn REQUIRED EXPORT )
find_package( larpandora REQUIRED EXPORT )
find_package( larpandoracontent REQUIRED EXPORT )
find_package( artdaq_core REQUIRED EXPORT )
find_package( sbnobj REQUIRED EXPORT )
find_package( sbnanaobj REQUIRED EXPORT )
find_package( sbndaq_artdaq_core REQUIRED EXPORT )
find_package( ifdh_art REQUIRED EXPORT )
find_package( ifbeam REQUIRED EXPORT )
find_package( libwda REQUIRED EXPORT )
find_package( dk2nudata REQUIRED EXPORT )
find_package( log4cpp REQUIRED EXPORT )
find_package( LibXml2 REQUIRED EXPORT )
find_package( CLHEP REQUIRED EXPORT )
find_package( ROOT REQUIRED EXPORT )
find_package( Geant4 REQUIRED EXPORT )
find_package( Boost COMPONENTS system filesystem REQUIRED EXPORT )
find_package( geant4reweight REQUIRED EXPORT )

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)
include(SBNutils)

add_subdirectory(sbncode)
add_subdirectory(fcl)
add_subdirectory(fclmodule)


#generate cmake config files
cet_cmake_config()
