#======================================================================
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

find_package(cetmodules 2.29.07 REQUIRED)
project(sbncode VERSION 09.60.01 LANGUAGES CXX)

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
find_package(art)
find_package(messagefacility)
find_package(art_root_io)
find_package(nusimdata)
find_package(nug4)
find_package(nugen)
find_package(nurandom)
find_package(nutools)
find_package(GENIE)
find_package(lardataobj)
find_package(larcoreobj)
find_package(larcorealg)
find_package(lardataalg)
find_package(larcore)
find_package(lardata)
find_package(larevt)
find_package(larsim)
find_package(larreco)
find_package(larpandora)
find_package(larpandoracontent)
find_package(artdaq_core )
find_package(sbnobj )
find_package(sbnanaobj )
find_package(sbndaq_artdaq_core )
find_package(sbndata)
find_package(ifdhc)
find_package(ifdh_art)
find_package(ifbeam)
find_package( trace )
find_package( libwda )

find_package(CLHEP REQUIRED)
find_package(ROOT REQUIRED)
find_package(Geant4 REQUIRED)
find_package(Boost COMPONENTS system filesystem REQUIRED)
find_package( GSL )
find_library( GSL NAMES gsl )
find_library( LOG4CPP NAMES log4cpp )
find_library( XML2 NAMES xml2 )


# GENIE libraries are needed, added Aug2021 by Keng
# genie (copied from nugen)
  cet_find_library( LOG4CPP      NAMES log4cpp      PATHS ENV LOG4CPP_LIB NO_DEFAULT_PATH )
  set( GFWMSG       ${GFwMsg_LIBRARY}        )
  set( GFWREG       ${GFwReg_LIBRARY}       )
  set( GFWALG       ${GFwAlg_LIBRARY}       )
  set( GFWINT       ${GFwInt_LIBRARY}        )
  set( GFWGHEP      ${GFwGHEP_LIBRARY}       )
  set( GFWNUM       ${GFwNum_LIBRARY}        )
  set( GFWUTL       ${GFwUtl_LIBRARY}        )
  set( GFWPARDAT    ${GFwParDat_LIBRARY}     )
  set( GFWEG        ${GFwEG_LIBRARY}         )
  set( GFWNTP       ${GFwNtp_LIBRARY}        )
  set( GPHXSIG      ${GPhXSIg_LIBRARY}       )
  set( GPHPDF       ${GPhPDF_LIBRARY}        )
  set( GPHNUCLST    ${GPhNuclSt_LIBRARY}     )
  set( GPHCMN       ${GPhCmn_LIBRARY}        )
  set( GPHDCY       ${GPhDcy_LIBRARY}        )
  set( GPHHADTRANSP ${GPhHadTransp_LIBRARY}  )
  set( GPHHADNZ     ${GPhHadnz_LIBRARY}      )
  set( GPHDEEX      ${GPhDeEx_LIBRARY}       )
  set( GPHAMNGXS    ${GPhAMNGXS_LIBRARY}     )
  set( GPHAMNGEG    ${GPhAMNGEG_LIBRARY}     )
  set( GPHCHMXS     ${GPhChmXS_LIBRARY}      )
  set( GPHCOHXS     ${GPhCohXS_LIBRARY}      )
  set( GPHCOHEG     ${GPhCohEG_LIBRARY}      )
  set( GPHDISXS     ${GPhDISXS_LIBRARY}      )
  set( GPHDISEG     ${GPhDISEG_LIBRARY}      )
  set( GPHDFRCXS    ${GPhDfrcXS_LIBRARY}     )
  set( GPHDFRCEG    ${GPhDfrcEG_LIBRARY}     )
  set( GPHGLWRESXS  ${GPhGlwResXS_LIBRARY}   )
  set( GPHGLWRESEG  ${GPhGlwResEG_LIBRARY}   )
  set( GPHIBDXS     ${GPhIBDXS_LIBRARY}      )
  set( GPHIBDEG     ${GPhIBDEG_LIBRARY}      )
  set( GPHMNUCXS    ${GPhMNucXS_LIBRARY}     )
  set( GPHMNUCEG    ${GPhMNucEG_LIBRARY}     )
  set( GPHMEL       ${GPhMEL_LIBRARY}        )
  set( GPHNUELXS    ${GPhNuElXS_LIBRARY}     )
  set( GPHNUELEG    ${GPhNuElEG_LIBRARY}     )
  set( GPHQELXS     ${GPhQELXS_LIBRARY}      )
  set( GPHQELEG     ${GPhQELEG_LIBRARY}      )
  set( GPHRESXS     ${GPhResXS_LIBRARY}      )
  set( GPHRESEG     ${GPhResEG_LIBRARY}      )
  set( GPHSTRXS     ${GPhStrXS_LIBRARY}      )
  set( GPHSTREG     ${GPhStrEG_LIBRARY}      )
  set( GPHNDCY      ${GPhNDcy_LIBRARY}       )
  set( GTLGEO       ${GTlGeo_LIBRARY}        )
  set( GTLFLX       ${GTlFlx_LIBRARY}        )
  set( GRWFWK       ${GRwFwk_LIBRARY}        )
  set( GRWIO        ${GRwIO_LIBRARY}         )
  set( GRWCLC       ${GRwClc_LIBRARY}        )

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
add_subdirectory(fcl)
add_subdirectory(fclmodule)

# packaging utility
cet_cmake_config()
