art_make(
  LIBRARY_NAME sbncode_SBNEventWeight
  LIB_LIBRARIES
    sbnobj_Common_SBNEventWeight
    sbncode_SBNEventWeight_Base
    sbncode_SBNEventWeight_Calculators_CrossSection
    sbncode_SBNEventWeight_Calculators_BNBFlux
    #sbncode_SBNEventWeight_Calculators_Geant4
  MODULE_LIBRARIES
    sbnobj_Common_SBNEventWeight
    sbncode_SBNEventWeight_Base
    sbncode_SBNEventWeight_Calculators_CrossSection
    sbncode_SBNEventWeight_Calculators_BNBFlux
    nusimdata_SimulationBase
    ${MF_MESSAGELOGGER}
    ${ROOT_GEOM}
    ${LXML2}
    GSL
    systematicstools_interface
    systematicstools_interpreters
    systematicstools_utility
)

install_headers()
install_fhicl()
install_source()

