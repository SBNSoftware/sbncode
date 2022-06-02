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

project(sbncode VERSION 09.53.00 LANGUAGES CXX)

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
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs
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
find_ups_product(larcoreobj)
find_ups_product(larcorealg)
find_ups_product(lardataobj)
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
find_ups_product(libwda) 
find_ups_product(sbnobj v09_10_00)
find_ups_product(sbndata v01_03)
find_ups_product(ifdhc)
find_ups_product(ifdh_art)
find_ups_product(log4cpp)
find_ups_boost()
find_ups_root()
#find_ups_product(sbndcode v06_67_00)
#find_ups_product(uboonecode v06_67_00)
find_ups_product( sbndaq_artdaq_core )


# GENIE libraries are needed, added Aug2021 by Keng 
# genie (copied from nugen)
  cet_find_library( LOG4CPP      NAMES log4cpp      PATHS ENV LOG4CPP_LIB NO_DEFAULT_PATH )
  cet_find_library( GFWMSG       NAMES GFwMsg       PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GFWREG       NAMES GFwReg       PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GFWALG       NAMES GFwAlg       PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GFWINT       NAMES GFwInt       PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GFWGHEP      NAMES GFwGHEP      PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GFWNUM       NAMES GFwNum       PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GFWUTL       NAMES GFwUtl       PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GFWPARDAT    NAMES GFwParDat    PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GFWEG        NAMES GFwEG        PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GFWNTP       NAMES GFwNtp       PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHXSIG      NAMES GPhXSIg      PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHPDF       NAMES GPhPDF       PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHNUCLST    NAMES GPhNuclSt    PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHCMN       NAMES GPhCmn       PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHDCY       NAMES GPhDcy       PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHHADTRANSP NAMES GPhHadTransp PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHHADNZ     NAMES GPhHadnz     PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHDEEX      NAMES GPhDeEx      PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHAMNGXS    NAMES GPhAMNGXS    PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHAMNGEG    NAMES GPhAMNGEG    PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHCHMXS     NAMES GPhChmXS     PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHCOHXS     NAMES GPhCohXS     PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHCOHEG     NAMES GPhCohEG     PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHDISXS     NAMES GPhDISXS     PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHDISEG     NAMES GPhDISEG     PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHDFRCXS    NAMES GPhDfrcXS    PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHDFRCEG    NAMES GPhDfrcEG    PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHGLWRESXS  NAMES GPhGlwResXS  PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHGLWRESEG  NAMES GPhGlwResEG  PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHIBDXS     NAMES GPhIBDXS     PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHIBDEG     NAMES GPhIBDEG     PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHMNUCXS    NAMES GPhMNucXS    PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHMNUCEG    NAMES GPhMNucEG    PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHMEL       NAMES GPhMEL       PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHNUELXS    NAMES GPhNuElXS    PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHNUELEG    NAMES GPhNuElEG    PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHQELXS     NAMES GPhQELXS     PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHQELEG     NAMES GPhQELEG     PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHRESXS     NAMES GPhResXS     PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHRESEG     NAMES GPhResEG     PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHSTRXS     NAMES GPhStrXS     PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHSTREG     NAMES GPhStrEG     PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GPHNDCY      NAMES GPhNDcy      PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GTLGEO       NAMES GTlGeo       PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GTLFLX       NAMES GTlFlx       PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GRWFWK       NAMES GRwFwk       PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GRWIO        NAMES GRwIO        PATHS ENV GENIE_LIB NO_DEFAULT_PATH )
  cet_find_library( GRWCLC       NAMES GRwClc       PATHS ENV GENIE_LIB NO_DEFAULT_PATH )

## (continued) GENIE pre-R-3 used GSL via ROOT's Math libs, now makes direct call not available via that route
cet_find_library( GSL NAMES gsl PATHS ENV GSL_LIB NO_DEFAULT_PATH )

  set (GENIE_LIB_LIST ${GFWMSG}
      ${GFWREG}
      ${GFWALG}
      ${GFWINT}
      ${GFWGHEP}
      ${GFWNUM}
      ${GSL}  # FWNUM relies on GSL
      ${GFWUTL}
      ${GFWPARDAT}
      ${GFWEG}
      ${GFWNTP}
      ${GPHXSIG}
      ${GPHPDF}
      ${GPHNUCLST}
      ${GPHCMN}
      ${GPHDCY}
      ${GPHHADTRANSP}
      ${GPHHADNZ}
      ${GPHDEEX}
      ${GPHAMNGXS}
      ${GPHAMNGEG}
      ${GPHCHMXS}
      ${GPHCOHXS}
      ${GPHCOHEG}
      ${GPHDISXS}
      ${GPHDISEG}
      ${GPHDFRCXS}
      ${GPHDFRCEG}
      ${GPHGLWRESXS}
      ${GPHGLWRESEG}
      ${GPHIBDXS}
      ${GPHIBDEG}
      ${GPHMNUCXS}
      ${GPHMNUCEG}
      ${GPHMEL}
      ${GPHNUELXS}
      ${GPHNUELEG}
      ${GPHQELXS}
      ${GPHQELEG}
      ${GPHRESXS}
      ${GPHRESEG}
      ${GPHSTRXS}
      ${GPHSTREG}
      ${GPHNDCY}
      ${GTLGEO}
      ${GTLFLX}
      ${GRWFWK}
      ${GRWIO}
      ${GRWCLC}
      ${LOG4CPP}
  )

# macros for dictionary and simple_plugin
include(ArtDictionary)
include(ArtMake)
include(BuildPlugins)

add_subdirectory(sbncode)
add_subdirectory(ups)
add_subdirectory(fcl)
add_subdirectory(fclmodule)

# packaging utility
include(UseCPack)
