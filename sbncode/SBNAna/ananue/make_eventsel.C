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
  const std::string fname = "/pnfs/icarus/persistent/users/dmendez/SBNAnaFiles/test/icarus/gen-prodcorsika_genie_nooverburden__nuetest.caf.root";

  // Source of events
  SpectrumLoader loader(fname);

  const unsigned int kNVar = plots.size();
  const unsigned int kNSel = sels.size();
  // const unsigned int kNVarSpill = plots_spill.size();
  // const unsigned int kNSelSpill = sels_spill.size();

  // Make an array of spectra with each Var-Cut combination
  Spectrum *specs[kNSel][kNVar];
  // Spectrum *specs_spill[kNSelSpill][kNVarSpill];

  for(unsigned int iSel = 0; iSel < kNSel; ++iSel){
    for(unsigned int jVar = 0; jVar < kNVar; ++jVar){
      specs[iSel][jVar] = new Spectrum(plots[jVar].label, plots[jVar].bins, loader, plots[jVar].var, sels[iSel].cut);
    }
  }
  // for(unsigned int iSel = 0; iSel < kNSelSpill; ++iSel){
  //   for(unsigned int jVar = 0; jVar < kNVarSpill; ++jVar){
  //     specs_spill[iSel][jVar] = new Spectrum(plots_spill[jVar].label, plots_spill[jVar].bins, loader, plots_spill[jVar].var, sels_spill[iSel].cut, sels[iSel].cut);
  //   }
  // }

  Spectrum *mixcut1 = new Spectrum("flashtrig", Binning::Simple(0,0,100), loader, kEvt, kFlashTrigger);
  Spectrum *mixcut2 = new Spectrum("noflashtrig", Binning::Simple(0,0,100), loader, kEvt, !kFlashTrigger);
  Spectrum *mixcut3 = new Spectrum("firstevts", Binning::Simple(0,0,100), loader, kEvt, kFirstEvents);
  Spectrum *mixcut4 = new Spectrum("nofirstevts", Binning::Simple(0,0,100), loader, kEvt, !kFirstEvents);
  // Spectrum *mixcut1 = new Spectrum("flashtrig_shortshw", kLengthBinning, loader, kRecoShower_Length, kFlashTrigger, kShortShower);
  // Spectrum *mixcut2 = new Spectrum("flashtrig_longshw", kLengthBinning, loader, kRecoShower_Length, kFlashTrigger, !kShortShower);
  // Spectrum *mixcut3 = new Spectrum("firstevt_shortshw", kLengthBinning, loader, kRecoShower_Length, kFirstEvents, kShortShower);
  // Spectrum *mixcut4 = new Spectrum("firstevt_longshw", kLengthBinning, loader, kRecoShower_Length, kFirstEvents, !kShortShower);

  // This is the call that actually fills in the spectra
  loader.Go();

  // Save spectra to a file, so we can plot them later
  TFile fout("spectra_eventsel_nue.root", "RECREATE");

  for( unsigned int iSel = 0; iSel < kNSel; ++iSel ){
    for( unsigned int jVar = 0; jVar < kNVar; ++jVar ){
      std::string mysuffix = sels[iSel].suffix + "_" + plots[jVar].suffix;
      std::cout << "Saving spectra: " << mysuffix << std::endl;
      specs[iSel][jVar]->SaveTo( fout.mkdir(mysuffix.c_str()) );
    }
  }

  // for( unsigned int iSel = 0; iSel < kNSelSpill; ++iSel ){
  //   for( unsigned int jVar = 0; jVar < kNVarSpill; ++jVar ){
  //     std::string mysuffix = sels_spill[iSel].suffix + "_" + plots_spill[jVar].suffix;
  //     std::cout << "Saving spectra: " << mysuffix << std::endl;
  //     specs_spill[iSel][jVar]->SaveTo( fout.mkdir(mysuffix.c_str()) );
  //   }
  // }

  mixcut1->SaveTo( fout.mkdir("mixcut1"));
  mixcut2->SaveTo( fout.mkdir("mixcut2"));
  mixcut3->SaveTo( fout.mkdir("mixcut3"));
  mixcut4->SaveTo( fout.mkdir("mixcut4"));

}
