# Headers needed for the GenieWeightCalc class
include_directories ( $ENV{GENIE_INC} )
include_directories ( $ENV{GENIE_INC}/GENIE )
include_directories ( $ENV{LOG4CPP_INC} )

#CHECK, how to make GENIE compiled?
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
    ${GENIE_LIB_LIST}
    ${LOG4CPP}
    ${MF_MESSAGELOGGER}
    cetlib_except
)

install_headers()
install_fhicl()
install_source()

