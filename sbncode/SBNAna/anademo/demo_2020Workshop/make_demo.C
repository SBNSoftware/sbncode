// Make a few spectra with different cuts.

#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"

#include "demo.h"

#include "TFile.h"

using namespace ana;

void make_demo()
{
  //
  // Environment variables and wildcards work as arguments to
  // SpectrumLoader. As do SAM datasets.
  const std::string fname = "larout.caf.root";
    //"/sbnd/app/users/psihas/SBNCAF-dev/prodgenie_bnb_nu_cosmic_sbnd_GenieGen-20191122T115936_G4-20191122T200111_DetSim-20191122T205333_Reco-20191122T224902_0c973425-93ec-4a22-9e8f-d2fe0fd2d189_ReOp-20200324T001103.caf.root";  
  ///sbnd/app/users/psihas/ExampleFiles/CAFs/prodgenie_bnb_nu_cosmic_sbnd_20191128_reco.caf.root";

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
