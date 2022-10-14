simple_plugin(SBNEventWeight module
  sbnobj::Common_SBNEventWeight sbncode::SBNEventWeight_Base
  sbncode::SBNEventWeight_Calculators_CrossSection
  sbncode::SBNEventWeight_Calculators_BNBFlux nusimdata::SimulationBase
  ROOT::Geom)

install_headers()
install_fhicl()
install_source()

