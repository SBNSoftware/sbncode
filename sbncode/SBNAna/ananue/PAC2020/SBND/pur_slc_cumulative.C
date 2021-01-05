#include "make_spectra_pur_slc_cumulative.C"
#include "plot_spectra_pur_slc_cumulative.C"

void pur_slc_cumulative(std::string outName = "")
{

  make_spectra_pur_slc_cumulative("/sbnd/data/users/etyley/cafs/2020A/NuEOverlayHadd.caf.root", "NuEOverlay" + outName);
  plot_spectra_pur_slc_cumulative("NuEOverlay" + outName);

  make_spectra_pur_slc_cumulative("/sbnd/data/users/etyley/cafs/2020A/NuMuOverlayHadd.caf.root", "NuMuOverlay" + outName);
  plot_spectra_pur_slc_cumulative("NuMuOverlay" + outName);

  make_spectra_pur_slc_cumulative("/sbnd/data/users/etyley/cafs/2020A/IntimeCosmicsHadd.caf.root", "IntimeCosmics" + outName);
  plot_spectra_pur_slc_cumulative("IntimeCosmics" + outName);
}
