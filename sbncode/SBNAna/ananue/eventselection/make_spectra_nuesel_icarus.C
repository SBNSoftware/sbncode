// Make a few spectra with different cuts.

#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"

#include "helper_nuesel_icarus.h"

#include "TFile.h"

using namespace ana;

void make_spectra_nuesel_icarus(std::string finname = "nucosmics", int setno = 1)
{

  std::string settag = std::to_string(setno);

  std::string findir = "cafsdir/"+finname+"/files/";
  const std::string finsample = findir + "hadded_"+settag+"_"+finname+".caf.root";
  const std::string foutname = finname+"_spectra_hadded_"+settag+".root";

  // Source of events
  SpectrumLoader loader(finsample);

  ////////////////////////////////
  //// pseudo scale cosmics
  TFile f(finsample.c_str());
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

  const float nuPerSpill(1.f / 21.8);
  const float POTperSpill(5e12);
  const float cosmicPOT(totalGenEvt * POTperSpill / (1.f - nuPerSpill));

  std::cout << "Spills per cosmic Trigger: " << (cosmicPOT / numTriggers) / POTperSpill << std::endl;
  //// end pseudo scale cosmics
  ////////////////////////////////

  const unsigned int kNVar = plots.size();
  const unsigned int kNSel = sels.size();
  const unsigned int kNType = types.size();
  const unsigned int kNVarSpill = plots_spill.size();
  const unsigned int kNSelSpill = sels_spill.size();

  // make an array of spectra with each cut-type-var combination
  Spectrum *specs[kNSel][kNType][kNVar];
  Spectrum *specsveto[kNSel][kNType][kNVar];
  Spectrum *specsspill[kNSelSpill][kNVarSpill];

  for(unsigned int iSel = 0; iSel < kNSel; ++iSel){
    for(unsigned int iType = 0; iType < kNType; ++iType){
      for(unsigned int jVar = 0; jVar < kNVar; ++jVar){
        specs[iSel][iType][jVar] = new Spectrum(plots[jVar].label, plots[jVar].bins, loader, plots[jVar].var, sels[iSel].cut && types[iType].cut);
        specsveto[iSel][iType][jVar] = new Spectrum(plots[jVar].label, plots[jVar].bins, loader, plots[jVar].var, kCRTHitVetoFD, sels[iSel].cut && types[iType].cut); // Add CRT veto
      }
    }
  }

  for(unsigned int iSel = 0; iSel < kNSelSpill; ++iSel){
    for(unsigned int jVar = 0; jVar < kNVarSpill; ++jVar){
      specsspill[iSel][jVar] = new Spectrum(plots_spill[jVar].label, plots_spill[jVar].bins, loader, plots_spill[jVar].var, sels_spill[iSel].cut);
    }
  }
  // actually make the spectra
  loader.Go();

  // save spectra to a file, so we can plot them later
  TFile fout(foutname.c_str(), "RECREATE");

  for( unsigned int iSel = 0; iSel < kNSel; ++iSel ){
    for(unsigned int iType = 0; iType < kNType; ++iType){
      for( unsigned int jVar = 0; jVar < kNVar; ++jVar ){
        std::string mysuffix = types[iType].suffix+"_"+sels[iSel].suffix + "_" + plots[jVar].suffix;
        std::string mysuffixveto = types[iType].suffix+"_"+sels[iSel].suffix + "_veto_" + plots[jVar].suffix;
        std::cout << "Saving spectra: " << mysuffix << std::endl;
        Spectrum* spec = specs[iSel][iType][jVar];
        Spectrum* specveto = specsveto[iSel][iType][jVar];
        if(finname=="cosmics"){
          std::cout << "fake POT scale cosmics" << std::endl;
          const double thisPOT(spec->POT());
          if (thisPOT < std::numeric_limits<double>::epsilon()){
            spec->OverridePOT(cosmicPOT);
            specveto->OverridePOT(cosmicPOT);
          }
        } // fake cosmics pot
        spec->SaveTo(fout.mkdir(mysuffix.c_str()));
        specveto->SaveTo(fout.mkdir(mysuffixveto.c_str()));
      }
    }
  }

  for( unsigned int iSel = 0; iSel < kNSelSpill; ++iSel ){
    for( unsigned int jVar = 0; jVar < kNVarSpill; ++jVar ){
      std::string mysuffix = sels_spill[iSel].suffix + "_" + plots_spill[jVar].suffix;
      std::cout << "Saving spectra: " << mysuffix << std::endl;
      specsspill[iSel][jVar]->SaveTo( fout.mkdir(mysuffix.c_str()) );
    }
  }

}
