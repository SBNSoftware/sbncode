// Make a few spectra with different cuts.

#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"

#include "helper.h"

#include "TFile.h"

using namespace ana;

void make_spectra()
{
  //
  // Environment variables and wildcards work as arguments to
  // SpectrumLoader. As do SAM datasets.
  // const std::string fname = "/pnfs/icarus/persistent/users/dmendez/SBNAnaFiles/test/icarus/gen-prodcorsika_genie_nooverburden__nuetest.caf.root";
  // const std::string fname = "/pnfs/icarus/persistent/users/dmendez/SBNAnaFiles/test/icarus/gen-prodcorsika_genie_nooverburden_slccopy.flatcaf.root";
  const std::string fname = "/pnfs/icarus/persistent/users/dmendez/SBNAnaFiles/test/icarus/gen-prodcorsika_genie_nooverburden_slccopy.caf.root";
  // Source of events
  SpectrumLoader loader(fname);

  const unsigned int kNVar = plots.size();
  const unsigned int kNSel = sels.size();
  const unsigned int kNVarSpill = plots_spill.size();
  const unsigned int kNSelSpill = sels_spill.size();

  // Make an array of spectra with each Var-Cut combination
  Spectrum *specs[kNSel][kNVar]; // (Slice)Var with (Slice)Cut

  for(unsigned int iSel = 0; iSel < kNSel; ++iSel){
    for(unsigned int jVar = 0; jVar < kNVar; ++jVar){
      specs[iSel][jVar] = new Spectrum(plots[jVar].label, plots[jVar].bins, loader, plots[jVar].var, sels[iSel].cut);
    }
  }

  // Example spectra with both SpillCut and Cut
  Spectrum *s0 = new Spectrum("crthitx__flashtrig", kPositionXFDBinning, loader, kCRTHitX, kFlashTrigger); // SpillMultiVar with SpillCut
  Spectrum *s1 = new Spectrum("shwlen__firstevt_longshw", kLengthBinning, loader, kRecoShower_Length, kFirstEvents, !kShortShower); // (Slice)Var with SpillCut and (Slice)Cut
  Spectrum *s2 = new Spectrum("shwlen__flashtrig_longshw", kLengthBinning, loader, kRecoShower_Length, kFlashTrigger, !kShortShower); // (Slice)Var with SpillCut and (Slice)Cut
  Spectrum *s3 = new Spectrum("shwlen__flashtrig_shortshw", kLengthBinning, loader, kRecoShower_Length, kFlashTrigger, kShortShower); // (Slice)Var with SpillCut and (Slice)Cut

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

  s0->SaveTo( fout.mkdir("s0"));
  s1->SaveTo( fout.mkdir("s1"));
  s2->SaveTo( fout.mkdir("s2"));
  s3->SaveTo( fout.mkdir("s3"));

}
