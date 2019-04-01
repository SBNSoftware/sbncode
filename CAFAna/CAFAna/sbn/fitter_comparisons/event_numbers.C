// Get Numbers for various interaction types and flavours
// cafe event_numbers.C

#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"

// #include "CAFAna/Cuts/TruthCuts.h"

#include "StandardRecord/StandardRecord.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Analysis/Calcs.h"
#include "CAFAna/Core/MultiVar.h"
#include "OscLib/func/OscCalculatorSterile.h"

#include "StandardRecord/Proxy/SRProxy.h"

#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "TLatex.h"

using namespace ana;

// Put an "SBND simulation" tag in the corner                                                       
void Experiment(std::string expt)
{
  TLatex* prelim = new TLatex(.9, .95, (expt+" Simulation").c_str());
  prelim->SetTextColor(kGray+1);
  prelim->SetNDC();
  prelim->SetTextSize(2/30.);
  prelim->SetTextAlign(32);
  prelim->Draw();
}

void event_numbers(const std::string expt = "SBND")
{
  const std::string fDir = "/pnfs/sbnd/persistent/users/gputnam/numu_simulation_reweight/processed_2.a/";
  const std::string fnameBeam = fDir + "output_SBNOsc_NumuSelection_Modern_" + expt + ".root";

  // Source of events
  SpectrumLoader loaderBeam(fnameBeam);
  // SpectrumLoader loaderSwap(fnameSwap);

  // Vars and cuts

  // a struct to hold the different flavours and signs
  const unsigned int kNumFlavours = 2;
  const unsigned int kNumSigns    = 2;
  const unsigned int kNumCurrents = 2;

  struct Flavour
  {
    ana::Flavors::Flavors_t flav;
    string label;
  };
  Flavour flav[kNumFlavours] = {
    {ana::Flavors::kAllNuMu, "numu"},
    {ana::Flavors::kAllNuE,  "nue"}
  };
  struct Signtype
  {
    ana::Sign::Sign_t sign;
    string label;
  };
    Signtype sign[kNumSigns] = {
    {ana::Sign::kNu,     ""},
    {ana::Sign::kAntiNu, "anti"}
  };
  struct Currtype
  {
    ana::Current::Current_t curr;
    string label;
  };
    Currtype curr[kNumCurrents] = {
    {ana::Current::kCC, "CC"},
    {ana::Current::kNC, "NC"}
  };

    const Var kSmearedE([](const caf::SRProxy* sr)
                          {
          return sr->reco[0].reco_energy;
        });

    const Var kWeight([](const caf::SRProxy* sr)
                          {
          return sr->reco[0].weight;
        });

  const Binning binsEnergy = Binning::Simple(30, 0, 3);
  const HistAxis axEnergy("Fake reconstructed energy (GeV)", binsEnergy, kSmearedE);

  const Var kCC = SIMPLEVAR(truth[0].neutrino.iscc);
  const Cut kIsCC = kCC > 0;

  // Temporary while we don't have the final states
  const Var kMode = SIMPLEVAR(truth[0].neutrino.genie_intcode);

  const Cut kIsQE  = (kMode == 0);
  const Cut kIsRes = (kMode == 1);
  const Cut kIsDIS = (kMode == 2);
  const Cut kIsCoh = (kMode == 3);
  const Cut kIsMEC = (kMode == 10);

  const unsigned int kNumModes = 6;
  struct Selection
  {
    Cut cut;
    string label;
    int colour;
  };
  Selection cuts[kNumModes] = {
    {kIsQE,  "QE", kRed},
    {kIsRes, "Res", kGreen+2},
    {kIsDIS, "DIS", kBlue},
    {kIsCoh, "Coh", kViolet},
    {kIsMEC, "MEC", kOrange+3},
    {kNoCut, "All", kBlack}
  };

  // For now let's not bother with oscillations
  osc::NoOscillations* noosc = new osc::NoOscillations;

  // The main prediction
  PredictionNoExtrap* pred[kNumModes];

  for(unsigned int m = 0; m < kNumModes; m++){
    //pred[m] = new PredictionNoExtrap(loaderBeam, loaderSwap, kNullLoader,
    pred[m] = new PredictionNoExtrap(loaderBeam, kNullLoader, kNullLoader,
      axEnergy, cuts[m].cut, kNoShift, kWeight);
  }

  // GO!
  loaderBeam.Go();
  // loaderSwap.Go();

  // Fake POT: we need to sort this out in the files first
  const double pot = 6.6e20;

  std::cout << "     NuMu,  Nue,  anti-Numu,  anti-Nue,  NC" << std::endl;

  TH1* h[kNumModes][kNumCurrents][kNumSigns][kNumFlavours];

  for(unsigned int m = 0; m < kNumModes; m++){

    std::cout << cuts[m].label;

    for(unsigned int k = 0; k < kNumCurrents; k++){
      for(unsigned int j = k; j < kNumSigns; j++){ // loop over signs only for CC
        for(unsigned int i = k; i < kNumFlavours; i++){ // loop over flavours only for CC
          h[m][k][j][i] = pred[m]->PredictComponent(noosc, k == 0 ? flav[i].flav : Flavors::kAll, 
                                    curr[k].curr, k == 0 ? sign[j].sign : Sign::kBoth).ToTH1(pot);
          float nEvts = h[m][k][j][i]->Integral();
          std::cout << ", " << nEvts ;
        } // end loop over flavours
      } // end loop over signs
    } // end loop over currents
    std::cout << endl;
  } // end loop over modes

  TFile* fOutput = new TFile(("output/output_"+expt+".root").c_str(),"RECREATE");

  TCanvas* c1 = new TCanvas("c1","c1");

  for(unsigned int m = 0; m < kNumModes; m++){
    for(unsigned int k = 0; k < kNumCurrents; k++){ 
      for(unsigned int j = k; j < kNumSigns; j++){ // loop over signs only for CC
        for(unsigned int i = k; i < kNumFlavours; i++){ // loop over flavours only for CC

          string currLabel = curr[k].label;
          string flavLabel = (k == 0 ? flav[i].label : "all");
          string signLabel = (k == 0 ? sign[j].label : "");

          h[m][k][j][i]->SetLineColor(cuts[m].colour);
          h[m][k][j][i]->SetMarkerColor(cuts[m].colour);
          h[m][k][j][i]->Write((expt+"_"+cuts[m].label+"_"+currLabel+"_"+signLabel+flavLabel).c_str());   
          h[m][k][j][i]->SetName((cuts[m].label+"_"+currLabel+"_"+signLabel+flavLabel).c_str());      
          h[m][k][j][i]->Draw("hist");
          Experiment(expt);
          c1->SaveAs(("output/"+expt+"_"+cuts[m].label+"_"+currLabel+"_"+signLabel+flavLabel+".pdf").c_str());

        } // loop over flavours
      } // loop over signs
    } // loop over currents
  } // loop over modes

  fOutput->Close();

}
