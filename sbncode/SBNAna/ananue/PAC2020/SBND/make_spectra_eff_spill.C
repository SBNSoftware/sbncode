// Make a few spectra with different cuts.

#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/SpectrumLoader.h"

#include "helper_eff_spill.h"

#include "TFile.h"
#include "TTree.h"

#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

using namespace ana;

void make_spectra_eff_spill(const std::string fname = "etyley_caf_NuEOverlay_caf_20201116_040705_447133", std::string name = "nue")
{

  // Source of events
  SpectrumLoader loader(fname);

  // Hack in some POT couniting for cosmics, should be done more properly in the future
  TFile f(fname.c_str());
  TTree* myTree = (TTree*)f.Get("recTree");

  TTreeReader cafReader("recTree", &f);

  TTreeReaderValue<unsigned int> run(cafReader, "hdr.run");
  TTreeReaderValue<unsigned int> subrun(cafReader, "hdr.subrun");
  TTreeReaderValue<unsigned int> ngenevt(cafReader, "hdr.ngenevt");

  unsigned int totalGenEvt(0);

  std::map<std::pair<unsigned int, unsigned int>, unsigned int> runSubRunMap;
  while (cafReader.Next()) {
    std::pair<unsigned int, unsigned int> thisPair(*run, *subrun);
    if (runSubRunMap[thisPair] == 0) {
      ++runSubRunMap[thisPair];
      totalGenEvt += *ngenevt;
    }
  }

  const unsigned int numTriggers(myTree->GetEntries());

  std::cout << "Num Triggers: " << numTriggers << " from generated events: " << totalGenEvt << std::endl;

  const float nuPerSpill(1.f / 21.8); // Stolen from Gray
  const float POTperSpill(5e12);
  const float cosmicPOT(totalGenEvt * POTperSpill / (1.f - nuPerSpill));

  std::cout << "Spills per cosmic Trigger: " << (cosmicPOT / numTriggers) / POTperSpill << std::endl;

  const unsigned int kNVar = plots.size();
  const unsigned int kNType = types.size();
  const unsigned int kNSel = sels.size();

  // Make an array of spectra with each Var-Cut combination
  Spectrum* specs[kNSel][kNType][kNVar]; // (Slice)Var with (Slice)Cut

  for (unsigned int iSel = 0; iSel < kNSel; ++iSel) {
    for (unsigned int jType = 0; jType < kNType; ++jType) {
      for (unsigned int lVar = 0; lVar < kNVar; ++lVar) {
        specs[iSel][jType][lVar] = new Spectrum(plots[lVar].label, plots[lVar].bins, loader, plots[lVar].var, types[jType].cut && sels[iSel].cut);
      }
    }
  }

  // This is the call that actually fills in the spectra
  loader.Go();

  // Save spectra to a file, so we can plot them later
  TFile fout(("spectra_eff_spill_" + name + ".root").c_str(), "RECREATE");

  for (unsigned int iSel = 0; iSel < kNSel; ++iSel) {
    for (unsigned int jType = 0; jType < kNType; ++jType) {
      for (unsigned int lVar = 0; lVar < kNVar; ++lVar) {
        std::string mysuffix = plots[lVar].suffix + "_" + sels[iSel].suffix + "_" + types[jType].suffix;
        std::cout << "Saving spectra: " << mysuffix << std::endl;
        Spectrum* spec = specs[iSel][jType][lVar];
        const double thisPOT(spec->POT());
        if (thisPOT < std::numeric_limits<double>::epsilon())
          spec->OverridePOT(cosmicPOT);
        spec->SaveTo(fout.mkdir(mysuffix.c_str()));
      }
    }
  }

  // for( unsigned int iSel = 0; iSel < kNSelSpill; ++iSel ){
  //   for( unsigned int jVar = 0; jVar < kNVarSpill; ++jVar ){
  //     std::string mysuffix = sels_spill[iSel].suffix + "_" + plots_spill[jVar].suffix;
  //     std::cout << "Saving spectra: " << mysuffix << std::endl;
  //     specsSpill[iSel][jVar]->SaveTo( fout.mkdir(mysuffix.c_str()) );
  //   }
  // }
}
