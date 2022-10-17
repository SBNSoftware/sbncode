# Headers needed for the GenieWeightCalc class
include_directories ( $ENV{GENIE_INC} )
include_directories ( $ENV{GENIE_INC}/GENIE )
include_directories ( $ENV{LOG4CPP_INC} )
include_directories($ENV{GENIE}/include/GENIE) 
include_directories( ${GSL_INCLUDE_DIR} )


art_make( 
  LIBRARY_NAME sbncode_SBNEventWeight_Calculators_CrossSection
  LIB_LIBRARIES
    sbncode::SBNEventWeight_Base
    nugen::EventGeneratorBase_GENIE
    nugen::NuReweight_art
    nugen::NuReweight
    nurandom::RandomUtils_NuRandomService_service
    art::Framework_Principal
    art::Framework_Services_Registry
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
    cetlib_except::cetlib_except
)

install_headers()
install_fhicl()
install_source()

