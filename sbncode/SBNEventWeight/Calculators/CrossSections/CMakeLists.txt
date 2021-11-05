# Headers needed for the GenieWeightCalc class
include_directories ( $ENV{GENIE_INC} )
include_directories ( $ENV{GENIE_INC}/GENIE )
include_directories ( $ENV{LOG4CPP_INC} )

# Builds of MicroBooNE-specific larsim releases should uncomment the
# "add_definitions" line below. It enables extra reweighting functionality
# which is present only in "MicroBooNE patch" releases of GENIE v3.0.4. These
# require building against a special ups product.

# Define a macro that indicates whether a GENIE version before 3.0.0 is in use
 parse_ups_version( ${GENIE_VERSION} )
 if( ${VMAJ} LESS 3 )
   add_definitions(-DGENIE_PRE_R3)
   endif()

# not applying the following, YET
#add_definitions(-DGENIE_UB_PATCH)

art_make(
  LIBRARY_NAME sbncode_SBNEventWeight_Calculators_CrossSection
  LIB_LIBRARIES
    sbncode_SBNEventWeight_Base
    nugen_EventGeneratorBase_GENIE
    nurandom_RandomUtils_NuRandomService_service
    ${ART_FRAMEWORK_PRINCIPAL}
    ${ART_PERSISTENCY_PROVENANCE}
    ${ART_FRAMEWORK_SERVICES_REGISTRY}
    ${FHICLCPP}

   #GENIE_LIB_LIST is explicitly listed below:
    LOG4CPP
    GFWMSG
    GFWREG
    GFWALG
    GFWINT
    GFWGHEP
    GFWNUM
    GFWUTL
    GFWPARDAT
    GFWEG
    GFWNTP
    GPHXSIG
    GPHPDF
    GPHNUCLST
    GPHCMN
    GPHDCY
    GPHHADTRANSP
    GPHHADNZ
    GPHDEEX
    GPHAMNGXS
    GPHAMNGEG
    GPHCHMXS
    GPHCOHXS
    GPHCOHEG
    GPHDISXS
    GPHDISEG
    GPHDFRCXS
    GPHDFRCEG
    GPHGLWRESXS
    GPHGLWRESEG
    GPHIBDXS
    GPHIBDEG
    GPHMNUCXS
    GPHMNUCEG
    GPHMEL
    GPHNUELXS
    GPHNUELEG
    GPHQELXS
    GPHQELEG
    GPHRESXS
    GPHRESEG
    GPHSTRXS
    GPHSTREG
    GPHNDCY
    GTLGEO
    GTLFLX
    GRWFWK
    GRWIO
    GRWCLC
	GSL
    #GENIE libraries are above
    ${MF_MESSAGELOGGER}
    cetlib_except
)

install_headers()
install_fhicl()
install_source()

