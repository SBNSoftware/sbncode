#include "make_spectra_eff_spill.C"
#include "plot_spectra_eff_spill.C"

void eff_spill(std::string outName = "")
{

  make_spectra_eff_spill("/sbnd/data/users/etyley/cafs/2020A/NuEOverlayHadd.caf.root", "NuEOverlay" + outName);
  plot_spectra_eff_spill("NuEOverlay" + outName);

  make_spectra_eff_spill("/sbnd/data/users/etyley/cafs/2020A/NuMuOverlayHadd.caf.root", "NuMuOverlay" + outName);
  plot_spectra_eff_spill("NuMuOverlay" + outName);
}
