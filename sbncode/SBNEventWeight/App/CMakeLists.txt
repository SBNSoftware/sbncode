cet_build_plugin(SBNEventWeight art::module
  LIBRARIES
  sbnobj::Common_SBNEventWeight
  sbncode_SBNEventWeight_Base
  sbncode_SBNEventWeight_Calculators_CrossSection
  sbncode_SBNEventWeight_Calculators_BNBFlux
  sbncode_SBNEventWeight_Calculators_Geant4
  nusimdata::SimulationBase
  ROOT::Geom)

cet_build_plugin(SystToolsEventWeight art::module
  LIBRARIES
  sbnobj::Common_SBNEventWeight sbncode_SBNEventWeight_Base
  systematicstools::interface
  systematicstools::interpreters
  systematicstools::utility)

install_headers()
install_fhicl()
install_source()

