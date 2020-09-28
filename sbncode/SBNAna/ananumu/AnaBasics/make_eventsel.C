// Make a few spectra with different cuts.

#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"

#include "eventsel.h"

#include "TFile.h"

using namespace ana;

void make_eventsel()
{
  //
  // Environment variables and wildcards work as arguments to
  // SpectrumLoader. As do SAM datasets.
  const std::string fname = "/sbnd/app/users/psihas/SBNCAF-dev/reco-a1e5cc95-e889-474e-87a6-ffc68ec3825f.caf.root";
  //  const std::string fname = "/pnfs/icarus/persistent/users/dmendez/SBNAnaFiles/workshopdemo/numu.caf.root";

  // Source of events
  SpectrumLoader loader(fname);

  const unsigned int kNVar = plots.size();
  const unsigned int kNSel = sels.size();

  // Make an array of spectra with each Var-Cut combination
  Spectrum *specs[kNSel][kNVar];

  for(unsigned int iSel = 0; iSel < kNSel; ++iSel){
    for(unsigned int jVar = 0; jVar < kNVar; ++jVar){
      specs[iSel][jVar] = new Spectrum(plots[jVar].label, plots[jVar].bins, loader, plots[jVar].var, sels[iSel].cut);
    }
  }

  // This is the call that actually fills in the spectra
  loader.Go();

  // Save spectra to a file, so we can plot them later
  TFile fout("out_demo.root", "RECREATE");

  for( unsigned int iSel = 0; iSel < kNSel; ++iSel ){
    for( unsigned int jVar = 0; jVar < kNVar; ++jVar ){
      std::string mysuffix = sels[iSel].suffix + "_" + plots[jVar].suffix;
      std::cout << "Saving spctra: " << mysuffix << std::endl;
      specs[iSel][jVar]->SaveTo( fout.mkdir(mysuffix.c_str()) );
    }
  }

}
