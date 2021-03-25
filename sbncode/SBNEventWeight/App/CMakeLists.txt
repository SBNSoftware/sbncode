art_make(
  LIBRARY_NAME sbncode_SBNEventWeight
  LIB_LIBRARIES
    sbnobj_Common_SBNEventWeight
    sbncode_SBNEventWeight_Base
    sbncode_SBNEventWeight_Calculators_CrossSection
    #sbncode_SBNEventWeight_Calculators_BNBFlux
    #sbncode_SBNEventWeight_Calculators_Geant4
  MODULE_LIBRARIES
    sbnobj_Common_SBNEventWeight
    sbncode_SBNEventWeight_Base
    sbncode_SBNEventWeight_Calculators_CrossSection
    nusimdata_SimulationBase
    ${MF_MESSAGELOGGER}
    ${ROOT_GEOM}
    ${LXML2}
    ${GSL}
)

install_headers()
install_fhicl()
install_source()

