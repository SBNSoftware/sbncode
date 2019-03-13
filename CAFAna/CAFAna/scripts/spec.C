// ETW May 2018
// Standard script for DUNE spectra

#include "CAFAna/Analysis/Fit.h"
#include "CAFAna/Analysis/CalcsNuFit.h"
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/LoadFromFile.h"
#include "CAFAna/Core/Loaders.h"
#include "CAFAna/Core/Progress.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Cuts/AnaCuts.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Prediction/PredictionInterp.h"
#include "CAFAna/Analysis/Exposures.h"
#include "CAFAna/Systs/Systs.h"

using namespace ana;

#include "Utilities/rootlogon.C"

#include "OscLib/func/IOscCalculator.h"

#include "StandardRecord/StandardRecord.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2.h"

const Var kRecoE_nue = SIMPLEVAR(dune.Ev_reco_nue);
const Var kRecoE_numu = SIMPLEVAR(dune.Ev_reco_numu);

// 125 MeV bins from 0.0 to 8GeV
const HistAxis axis_nue("Reconstructed energy (GeV)",
                    Binning::Simple(64, 0.0, 8.0),
                    kRecoE_nue);

const HistAxis axis_numu("Reconstructed energy (GeV)",
                    Binning::Simple(64, 0.0, 8.0),
                    kRecoE_numu);


// POT/yr * 3.5yrs * mass correction
const double potFD = 3.5 * POT120 * 40/1.13;

const char* stateFname = "spec_state.root";
const char* outputFname = "spec_hist.root";

//Set systematics style by hand for now
bool nosyst = false;
bool normsyst = true;
bool fullsyst = false;

std::vector<const ISyst*> systlist;
std::vector<const ISyst*> normlist_sig = {&kNueFHCSyst, &kNumuFHCSyst, &kNueRHCSyst, &kNumuRHCSyst};
std::vector<const ISyst*> normlist_bg = {&kNueBeamFHCSyst, &kNueBeamRHCSyst, &kNCAppSyst, &kNCDisSyst, &kNutauSyst};

void spec(bool reload = false)
{
  rootlogon(); // style
  
  if(reload || TFile(stateFname).IsZombie()){

    SpectrumLoader loaderFDFHCNonswap("/dune/data/users/marshalc/CAFs/mcc11_v1/FD_FHC_nonswap.root");
    SpectrumLoader loaderFDFHCNue("/dune/data/users/marshalc/CAFs/mcc11_v1/FD_FHC_nueswap.root");
    SpectrumLoader loaderFDFHCNuTau("/dune/data/users/marshalc/CAFs/mcc11_v1/FD_FHC_tauswap.root");

    SpectrumLoader loaderFDRHCNonswap("/dune/data/users/marshalc/CAFs/mcc11_v1/FD_RHC_nonswap.root");
    SpectrumLoader loaderFDRHCNue("/dune/data/users/marshalc/CAFs/mcc11_v1/FD_RHC_nueswap.root");
    SpectrumLoader loaderFDRHCNuTau("/dune/data/users/marshalc/CAFs/mcc11_v1/FD_RHC_tauswap.root");

    Loaders dummyLoaders; // PredictionGenerator insists on this
    osc::IOscCalculatorAdjustable* calc = NuFitOscCalc(1);
    //Note that systlist must be filled both here and after the state load
    if (normsyst) {
      systlist.insert(systlist.end(), normlist_sig.begin(), normlist_sig.end()); 
      systlist.insert(systlist.end(), normlist_bg.begin(), normlist_bg.end()); 
    }
    // List all of the systematics we'll be using
    std::cout << "Systematics in these spectra: " << std::endl;
    for(const ISyst* s: systlist) std::cout << s->ShortName() << "\t\t" << s->LatexName() << std::endl;
    if (systlist.size()==0) {std::cout << "None" << std::endl;}


    NoExtrapPredictionGenerator genFDNumu(axis_numu, kPassFD_CVN_NUMU && kIsTrueFV);
    NoExtrapPredictionGenerator genFDNue(axis_nue, kPassFD_CVN_NUE && kIsTrueFV);
    NoExtrapPredictionGenerator genFDNumu_fid(axis_numu, kIsTrueFV);
    NoExtrapPredictionGenerator genFDNue_fid(axis_nue, kIsTrueFV);
    
    Loaders FD_FHC_loaders;
    Loaders FD_RHC_loaders;
    FD_FHC_loaders .AddLoader(&loaderFDFHCNonswap,  caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNonSwap);
    FD_FHC_loaders .AddLoader(&loaderFDFHCNue,      caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNueSwap);
    FD_FHC_loaders .AddLoader(&loaderFDFHCNuTau,    caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNuTauSwap);
    FD_RHC_loaders .AddLoader(&loaderFDRHCNonswap,  caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNonSwap);
    FD_RHC_loaders .AddLoader(&loaderFDRHCNue,      caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNueSwap);
    FD_RHC_loaders .AddLoader(&loaderFDRHCNuTau,    caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNuTauSwap);
      
    PredictionInterp predInt_FDNumuFHC(systlist, calc, genFDNumu, FD_FHC_loaders);
    PredictionInterp predInt_FDNumuRHC(systlist, calc, genFDNumu, FD_RHC_loaders);
    PredictionInterp predInt_FDNueFHC(systlist, calc, genFDNue, FD_FHC_loaders);
    PredictionInterp predInt_FDNueRHC(systlist, calc, genFDNue, FD_RHC_loaders);
    PredictionInterp predInt_FDNumuFHC_fid(systlist, calc, genFDNumu_fid, FD_FHC_loaders);
    PredictionInterp predInt_FDNumuRHC_fid(systlist, calc, genFDNumu_fid, FD_RHC_loaders);
    PredictionInterp predInt_FDNueFHC_fid(systlist, calc, genFDNue_fid, FD_FHC_loaders);
    PredictionInterp predInt_FDNueRHC_fid(systlist, calc, genFDNue_fid, FD_RHC_loaders);

    loaderFDFHCNonswap.Go();
    loaderFDFHCNue.Go();
    loaderFDFHCNuTau.Go();
    loaderFDRHCNonswap.Go();
    loaderFDRHCNue.Go();
    loaderFDRHCNuTau.Go();

    std::cout << stateFname << std::endl;
    TFile fout(stateFname, "RECREATE");
    predInt_FDNumuFHC.SaveTo(fout.mkdir("pred_fd_numu_fhc"));
    predInt_FDNumuRHC.SaveTo(fout.mkdir("pred_fd_numu_rhc"));
    predInt_FDNueFHC.SaveTo(fout.mkdir("pred_fd_nue_fhc"));
    predInt_FDNueRHC.SaveTo(fout.mkdir("pred_fd_nue_rhc"));
    predInt_FDNumuFHC_fid.SaveTo(fout.mkdir("pred_fd_numu_fhc_fid"));
    predInt_FDNumuRHC_fid.SaveTo(fout.mkdir("pred_fd_numu_rhc_fid"));
    predInt_FDNueFHC_fid.SaveTo(fout.mkdir("pred_fd_nue_fhc_fid"));
    predInt_FDNueRHC_fid.SaveTo(fout.mkdir("pred_fd_nue_rhc_fid"));

    std::cout << "All done making state..." << std::endl;

    }
  else{    
    std::cout << "Loading state from " << stateFname << std::endl; 
  }

  TFile fin(stateFname);
  PredictionInterp& predInt_FDNumuFHC = *ana::LoadFrom<PredictionInterp>(fin.GetDirectory("pred_fd_numu_fhc")).release();
  PredictionInterp& predInt_FDNueFHC = *ana::LoadFrom<PredictionInterp>(fin.GetDirectory("pred_fd_nue_fhc")).release();
  PredictionInterp& predInt_FDNumuRHC = *ana::LoadFrom<PredictionInterp>(fin.GetDirectory("pred_fd_numu_rhc")).release();
  PredictionInterp& predInt_FDNueRHC = *ana::LoadFrom<PredictionInterp>(fin.GetDirectory("pred_fd_nue_rhc")).release();
  PredictionInterp& predInt_FDNumuFHC_fid = *ana::LoadFrom<PredictionInterp>(fin.GetDirectory("pred_fd_numu_fhc_fid")).release();
  PredictionInterp& predInt_FDNueFHC_fid  = *ana::LoadFrom<PredictionInterp>(fin.GetDirectory("pred_fd_nue_fhc_fid")).release();
  PredictionInterp& predInt_FDNumuRHC_fid = *ana::LoadFrom<PredictionInterp>(fin.GetDirectory("pred_fd_numu_rhc_fid")).release();
  PredictionInterp& predInt_FDNueRHC_fid = *ana::LoadFrom<PredictionInterp>(fin.GetDirectory("pred_fd_nue_rhc_fid")).release();

  fin.Close();
  std::cout << "Done loading state" << std::endl;

  TFile* fout = new TFile(outputFname, "RECREATE");

  osc::NoOscillations noOsc;
  TH1* hnumufhc_sig_unosc = predInt_FDNumuFHC.PredictComponent(&noOsc, Flavors::kAllNuMu, Current::kCC, Sign::kNu).ToTH1(potFD);
  TH1* hnumufhc_sig_unosc_fid = predInt_FDNumuFHC_fid.PredictComponent(&noOsc, Flavors::kAllNuMu, Current::kCC, Sign::kNu).ToTH1(potFD);

  std::string dcpnames[4] = {"0pi","piover2","pi","3piover2"};

  for(int hie = -1; hie <= +1; hie += 2){

    osc::IOscCalculatorAdjustable* inputOsc = NuFitOscCalc(hie);

    const std::string hieStr = (hie > 0) ? "nh" : "ih";

    for(int deltaIdx = 0; deltaIdx < 4; ++deltaIdx){
      inputOsc->SetdCP(deltaIdx/2.*TMath::Pi());
      const std::string dcpStr = dcpnames[deltaIdx];

      TH1* hnumufhc = predInt_FDNumuFHC.Predict(inputOsc).ToTH1(potFD);
      TH1* hnuefhc = predInt_FDNueFHC.Predict(inputOsc).ToTH1(potFD);
      TH1* hnumurhc = predInt_FDNumuRHC.Predict(inputOsc).ToTH1(potFD);
      TH1* hnuerhc = predInt_FDNueRHC.Predict(inputOsc).ToTH1(potFD);

      //FHC Dis
      TH1* hnumufhc_sig = predInt_FDNumuFHC.PredictComponent(inputOsc, Flavors::kAllNuMu, Current::kCC, Sign::kNu).ToTH1(potFD);
      TH1* hnumufhc_ncbg = predInt_FDNumuFHC.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(potFD);
      TH1* hnumufhc_nutaubg = predInt_FDNumuFHC.PredictComponent(inputOsc, Flavors::kAllNuTau, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnumufhc_wsbg = predInt_FDNumuFHC.PredictComponent(inputOsc, Flavors::kAllNuMu, Current::kCC, Sign::kAntiNu).ToTH1(potFD);
      TH1* hnumufhc_nuebg = predInt_FDNumuFHC.PredictComponent(inputOsc, Flavors::kAllNuE, Current::kCC, Sign::kBoth).ToTH1(potFD);

      
      TH1* hnumufhc_sig_fid = predInt_FDNumuFHC_fid.PredictComponent(inputOsc, Flavors::kAllNuMu, Current::kCC, Sign::kNu).ToTH1(potFD);
      TH1* hnumufhc_ncbg_fid = predInt_FDNumuFHC_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(potFD);
      TH1* hnumufhc_nutaubg_fid = predInt_FDNumuFHC_fid.PredictComponent(inputOsc, Flavors::kAllNuTau, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnumufhc_wsbg_fid = predInt_FDNumuFHC_fid.PredictComponent(inputOsc, Flavors::kAllNuMu, Current::kCC, Sign::kAntiNu).ToTH1(potFD);
      TH1* hnumufhc_nuebg_fid = predInt_FDNumuFHC_fid.PredictComponent(inputOsc, Flavors::kAllNuE, Current::kCC, Sign::kBoth).ToTH1(potFD);

      //RHC Dis
      TH1* hnumurhc_sig = predInt_FDNumuRHC.PredictComponent(inputOsc, Flavors::kAllNuMu, Current::kCC, Sign::kAntiNu).ToTH1(potFD);
      TH1* hnumurhc_ncbg = predInt_FDNumuRHC.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(potFD);
      TH1* hnumurhc_nutaubg = predInt_FDNumuRHC.PredictComponent(inputOsc, Flavors::kAllNuTau, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnumurhc_wsbg = predInt_FDNumuRHC.PredictComponent(inputOsc, Flavors::kAllNuMu, Current::kCC, Sign::kNu).ToTH1(potFD);

      TH1* hnumurhc_sig_fid = predInt_FDNumuRHC_fid.PredictComponent(inputOsc, Flavors::kAllNuMu, Current::kCC, Sign::kAntiNu).ToTH1(potFD);
      TH1* hnumurhc_ncbg_fid = predInt_FDNumuRHC_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(potFD);
      TH1* hnumurhc_nutaubg_fid = predInt_FDNumuRHC_fid.PredictComponent(inputOsc, Flavors::kAllNuTau, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnumurhc_wsbg_fid = predInt_FDNumuRHC_fid.PredictComponent(inputOsc, Flavors::kAllNuMu, Current::kCC, Sign::kNu).ToTH1(potFD);

      //FHC App
      TH1* hnuefhc_sig = predInt_FDNueFHC.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuefhc_ncbg = predInt_FDNueFHC.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuefhc_beambg = predInt_FDNueFHC.PredictComponent(inputOsc, Flavors::kNuEToNuE, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuefhc_nutaubg = predInt_FDNueFHC.PredictComponent(inputOsc, Flavors::kAllNuTau, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuefhc_numubg = predInt_FDNueFHC.PredictComponent(inputOsc, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(potFD);

      TH1* hnuefhc_sig_fid = predInt_FDNueFHC_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuefhc_ncbg_fid = predInt_FDNueFHC_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuefhc_beambg_fid = predInt_FDNueFHC_fid.PredictComponent(inputOsc, Flavors::kNuEToNuE, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuefhc_nutaubg_fid = predInt_FDNueFHC_fid.PredictComponent(inputOsc, Flavors::kAllNuTau, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuefhc_numubg_fid = predInt_FDNueFHC_fid.PredictComponent(inputOsc, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(potFD);

      //RHC App
      TH1* hnuerhc_sig = predInt_FDNueRHC.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuerhc_ncbg = predInt_FDNueRHC.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuerhc_beambg = predInt_FDNueRHC.PredictComponent(inputOsc, Flavors::kNuEToNuE, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuerhc_nutaubg = predInt_FDNueRHC.PredictComponent(inputOsc, Flavors::kAllNuTau, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuerhc_numubg = predInt_FDNueRHC.PredictComponent(inputOsc, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(potFD);

      TH1* hnuerhc_sig_fid = predInt_FDNueRHC_fid.PredictComponent(inputOsc, Flavors::kNuMuToNuE, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuerhc_ncbg_fid = predInt_FDNueRHC_fid.PredictComponent(inputOsc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuerhc_beambg_fid = predInt_FDNueRHC_fid.PredictComponent(inputOsc, Flavors::kNuEToNuE, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuerhc_nutaubg_fid = predInt_FDNueRHC_fid.PredictComponent(inputOsc, Flavors::kAllNuTau, Current::kCC, Sign::kBoth).ToTH1(potFD);
      TH1* hnuerhc_numubg_fid = predInt_FDNueRHC_fid.PredictComponent(inputOsc, Flavors::kAllNuMu, Current::kCC, Sign::kBoth).ToTH1(potFD);

      hnumufhc_sig_unosc->Write("numu_fhc_unosc_sig");
      hnumufhc_sig_unosc_fid->Write("fid_numu_fhc_unosc_sig");

      hnumufhc->Write(("numu_fhc_"+hieStr+"_"+dcpStr).c_str());
      hnumufhc_sig->Write(("numu_fhc_sig_"+hieStr+"_"+dcpStr).c_str());
      hnumufhc_ncbg->Write(("numu_fhc_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnumufhc_nutaubg->Write(("numu_fhc_nutaubg_"+hieStr+"_"+dcpStr).c_str());
      hnumufhc_nuebg->Write(("numu_fhc_nuebg_"+hieStr+"_"+dcpStr).c_str());
      hnumufhc_wsbg->Write(("numu_fhc_wsbg_"+hieStr+"_"+dcpStr).c_str());

      hnuefhc->Write(("nue_fhc_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_sig->Write(("nue_fhc_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_ncbg->Write(("nue_fhc_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_beambg->Write(("nue_fhc_beambg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_nutaubg->Write(("nue_fhc_nutaubg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_numubg->Write(("nue_fhc_numubg_"+hieStr+"_"+dcpStr).c_str());
	
      hnumurhc->Write(("numu_rhc_"+hieStr+"_"+dcpStr).c_str());
      hnumurhc_sig->Write(("numu_rhc_sig_"+hieStr+"_"+dcpStr).c_str());
      hnumurhc_ncbg->Write(("numu_rhc_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnumurhc_nutaubg->Write(("numu_rhc_nutaubg_"+hieStr+"_"+dcpStr).c_str());
      hnumurhc_wsbg->Write(("numu_rhc_wsbg_"+hieStr+"_"+dcpStr).c_str());

      hnuerhc->Write(("nue_rhc_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_sig->Write(("nue_rhc_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_ncbg->Write(("nue_rhc_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_beambg->Write(("nue_rhc_beambg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_nutaubg->Write(("nue_rhc_nutaubg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_numubg->Write(("nue_rhc_numubg_"+hieStr+"_"+dcpStr).c_str());

      hnumufhc_sig_fid->Write(("fid_numu_fhc_sig_"+hieStr+"_"+dcpStr).c_str());
      hnumufhc_ncbg_fid->Write(("fid_numu_fhc_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnumufhc_nutaubg_fid->Write(("fid_numu_fhc_nutaubg_"+hieStr+"_"+dcpStr).c_str());
      hnumufhc_nuebg_fid->Write(("fid_numu_fhc_nuebg_"+hieStr+"_"+dcpStr).c_str());
      hnumufhc_wsbg_fid->Write(("fid_numu_fhc_wsbg_"+hieStr+"_"+dcpStr).c_str());

      hnuefhc_sig_fid->Write(("fid_nue_fhc_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_ncbg_fid->Write(("fid_nue_fhc_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_beambg_fid->Write(("fid_nue_fhc_beambg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_nutaubg_fid->Write(("fid_nue_fhc_nutaubg_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc_numubg_fid->Write(("fid_nue_fhc_numubg_"+hieStr+"_"+dcpStr).c_str());
	
      hnumurhc_sig_fid->Write(("fid_numu_rhc_sig_"+hieStr+"_"+dcpStr).c_str());
      hnumurhc_ncbg_fid->Write(("fid_numu_rhc_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnumurhc_nutaubg_fid->Write(("fid_numu_rhc_nutaubg_"+hieStr+"_"+dcpStr).c_str());
      hnumurhc_wsbg_fid->Write(("fid_numu_rhc_wsbg_"+hieStr+"_"+dcpStr).c_str());

      hnuerhc_sig_fid->Write(("fid_nue_rhc_sig_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_ncbg_fid->Write(("fid_nue_rhc_ncbg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_beambg_fid->Write(("fid_nue_rhc_beambg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_nutaubg_fid->Write(("fid_nue_rhc_nutaubg_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc_numubg_fid->Write(("fid_nue_rhc_numubg_"+hieStr+"_"+dcpStr).c_str());

    }
  }
  fout->Close();
}

