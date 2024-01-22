art_make_library(
  LIBRARY_NAME sbncode_SBNEventWeight_Calculators_MultiVariateRNG
  LIB_LIBRARIES
    fhiclcpp::fhiclcpp
    CLHEP::Random
    ROOT::Core
    cetlib::cetlib
    cetlib_except::cetlib_except
    ROOT::Graf3d
)

install_headers()
install_source()
