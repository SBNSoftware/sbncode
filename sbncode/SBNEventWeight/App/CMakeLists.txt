cet_build_plugin(SBNEventWeight art::module
  LIBRARIES
  sbnobj::Common_SBNEventWeight sbncode::SBNEventWeight_Base
  sbncode::SBNEventWeight_Calculators_CrossSection
  sbncode::SBNEventWeight_Calculators_BNBFlux nusimdata::SimulationBase
  ROOT::Geom)

install_headers()
install_fhicl()
install_source()

