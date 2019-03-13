#include "CAFAna/Experiment/MultiExperiment.h"
#include "CAFAna/Experiment/SingleSampleExperiment.h"
#include "CAFAna/Prediction/PredictionInterp.h"
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Prediction/PredictionNoOsc.h"
#include "CAFAna/Analysis/Calcs.h"
#include "OscLib/func/IOscCalculator.h"
#include "OscLib/func/OscCalculatorPMNSOpt.h"
#include "StandardRecord/StandardRecord.h"
#include "TCanvas.h"
#include "CAFAna/Systs/Systs.h"
#include "CAFAna/Systs/DUNEFluxSysts.h"
#include "CAFAna/Systs/GenieSysts.h"
#include "CAFAna/Systs/EnergySysts.h"
#include "CAFAna/Systs/NDRecoSysts.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "CAFAna/Analysis/Plots.h"
#include "CAFAna/Vars/FitVars.h"

#include "CAFAna/Analysis/Fit.h"
#include "CAFAna/Analysis/CalcsNuFit.h"

#include "CAFAna/Core/LoadFromFile.h"
#include "CAFAna/Core/Loaders.h"
#include "CAFAna/Core/Progress.h"

#include "CAFAna/Analysis/Surface.h"
#include "CAFAna/Analysis/Exposures.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Cuts/AnaCuts.h"
#include <tuple>
#include "Utilities/rootlogon.C"

using namespace ana;

// List of vars
// -->FD
const Var kRecoE_nue  = SIMPLEVAR(dune.Ev_reco_nue);
const Var kRecoE_numu = SIMPLEVAR(dune.Ev_reco_numu);
const Var kFDNumuPid  = SIMPLEVAR(dune.cvnnumu);
const Var kFDNuePid   = SIMPLEVAR(dune.cvnnue);
const Var kMVANUMU    = SIMPLEVAR(dune.mvanumu);

// -->ND
const Var kRecoEnergyND = SIMPLEVAR(dune.Ev_reco);
const Var kRecoYND      = (SIMPLEVAR(dune.Ev_reco) - SIMPLEVAR(dune.Elep_reco))/SIMPLEVAR(dune.Ev_reco);

// CV weighting
const Var kGENIEWeights = SIMPLEVAR(dune.total_cv_wgt); // kUnweighted

// FD cut
const Cut kFDSelNue     = kFDNuePid > 0.7;
const Cut kFDSelNumu    = kFDNumuPid > 0.7;
  
// --> ND cuts, from Chris: For the numu sample: reco_numu ==1, reco_q == -1 for FHC and +1 for RHC.  Also muon_exit == 0, which means that the muon is well-reconstructed.  And Ehad_veto < 30, which means the hadronic system is (probably) well-reconstructed
const Cut kRecoNegMu    = SIMPLEVAR(dune.reco_q) == -1; // Note that for these to be true, reco_numu == 1
const Cut kRecoPosMu    = SIMPLEVAR(dune.reco_q) == +1; // reco_q == 0 if reco_numu != 1
const Cut kMuonCont     = SIMPLEVAR(dune.muon_exit) == 0;
const Cut kEhad_veto    = SIMPLEVAR(dune.Ehad_veto) < 30;

// Binnings
const Binning binsFDEreco = Binning::Simple(80, 0, 10);
const Binning binsNDEreco = Binning::Simple(40, 0, 10);
const Binning binsY       = Binning::Simple(5, 0, 1);
						      
// Axes
const HistAxis axRecoEnuFDnumu("Reco energy (GeV)", binsFDEreco, kRecoE_numu);
const HistAxis axRecoEnuFDnue("Reco energy (GeV)", binsFDEreco, kRecoE_nue);
const HistAxis axErecYrecND("Reco energy (GeV)", binsNDEreco, kRecoEnergyND, "y_{rec}", binsY, kRecoYND);

// POT for 3.5 years
const double pot_fd = 3.5 * POT120 * 40/1.13;
const double pot_nd = 3.5 * POT120;

// FD spectra
SpectrumLoader FD_loaderFHCNumu("/dune/data/users/marshalc/CAFs/mcc11_v2/FD_FHC_nonswap.root", kBeam); //,1e5);
SpectrumLoader FD_loaderFHCNue("/dune/data/users/marshalc/CAFs/mcc11_v2/FD_FHC_nueswap.root", kBeam); //,1e5);
SpectrumLoader FD_loaderFHCNutau("/dune/data/users/marshalc/CAFs/mcc11_v2/FD_FHC_tauswap.root", kBeam); //,1e5);

SpectrumLoader FD_loaderRHCNumu("/dune/data/users/marshalc/CAFs/mcc11_v2/FD_RHC_nonswap.root", kBeam); //,1e5);
SpectrumLoader FD_loaderRHCNue("/dune/data/users/marshalc/CAFs/mcc11_v2/FD_RHC_nueswap.root", kBeam); //,1e5);
SpectrumLoader FD_loaderRHCNutau("/dune/data/users/marshalc/CAFs/mcc11_v2/FD_RHC_tauswap.root", kBeam); //,1e5);

// FD predictions
PredictionNoExtrap FD_predFHCNumu(FD_loaderFHCNumu, FD_loaderFHCNue, FD_loaderFHCNutau, axRecoEnuFDnumu, kPassFD_CVN_NUMU && kIsTrueFV, kNoShift, kGENIEWeights);
PredictionNoExtrap FD_predFHCNue (FD_loaderFHCNumu, FD_loaderFHCNue, FD_loaderFHCNutau, axRecoEnuFDnue , kPassFD_CVN_NUE && kIsTrueFV, kNoShift, kGENIEWeights);
PredictionNoExtrap FD_predRHCNumu(FD_loaderRHCNumu, FD_loaderRHCNue, FD_loaderRHCNutau, axRecoEnuFDnumu, kPassFD_CVN_NUMU && kIsTrueFV, kNoShift, kGENIEWeights);
PredictionNoExtrap FD_predRHCNue (FD_loaderRHCNumu, FD_loaderRHCNue, FD_loaderRHCNutau, axRecoEnuFDnue , kPassFD_CVN_NUE && kIsTrueFV, kNoShift, kGENIEWeights);

// ND predictions
SpectrumLoader ND_loaderFHC("/dune/data/users/marshalc/CAFs/mcc11_v2/ND_FHC_CAF.root", kBeam); //,1e5);
SpectrumLoader ND_loaderRHC("/dune/data/users/marshalc/CAFs/mcc11_v2/ND_RHC_CAF.root", kBeam); //,1e5);
PredictionNoOsc ND_predFHC(ND_loaderFHC, axErecYrecND, kPassND_FHC_NUMU && kIsTrueFV, kNoShift, kGENIEWeights);
PredictionNoOsc ND_predRHC(ND_loaderRHC, axErecYrecND, kPassND_RHC_NUMU && kIsTrueFV, kNoShift, kGENIEWeights);
// PredictionNoOsc ND_predFHC(ND_loaderFHC, axErecNPions, kPassND_FHC_NUMU && kIsTrueFV, kNoShift, kGENIEWeights);
// PredictionNoOsc ND_predRHC(ND_loaderRHC, axErecNPions, kPassND_RHC_NUMU && kIsTrueFV, kNoShift, kGENIEWeights);

// For the ND prediction generator
Loaders dummyLoaders;

// To get the oscillation probabilities
osc::IOscCalculatorAdjustable* calc = DefaultOscCalc();

std::vector<const ISyst*> detlist_nd  = {&kEnergyScaleMuSystND, &kChargedHadUncorrNDSyst, &kNUncorrNDSyst,
					 &kPi0UncorrNDSyst, &kRecoNCSyst, &kLeptonAccSyst, &kHadronAccSyst};
std::vector<const ISyst*> detlist_fd  = {&keScaleMuLArSyst, &kEnergyScaleESyst,
					 &kChargedHadCorrSyst, &kChargedHadUncorrFDSyst,
					 &kNUncorrFDSyst, &kEnergyScalePi0Syst,
					 &kPi0UncorrFDSyst};

std::vector<const ISyst*> detlist_dis  = {&keScaleMuLArSyst, &kChargedHadCorrSyst,
					  &kChargedHadUncorrFDSyst, &kNUncorrFDSyst,
					  &kEnergyScalePi0Syst, &kPi0UncorrFDSyst};
std::vector<const ISyst*> detlist_app  = {&kEnergyScaleESyst, &kChargedHadCorrSyst,
					  &kChargedHadUncorrFDSyst, &kNUncorrFDSyst,
					  &kEnergyScalePi0Syst, &kPi0UncorrFDSyst};

enum DetSystType {kFDAll, kFDApp, kFDDis};

// This is getting a little complicated
std::vector<const ISyst*> GetListOfSysts(bool fluxsyst, bool xsecsyst, bool detsyst,
					 bool fluxXsecPenalties = true, bool useND=true,
					 DetSystType det=kFDAll){

  std::vector<const ISyst*> systlist;
  if (fluxsyst){
    std::vector<const ISyst*> fluxlist = GetDUNEFluxSysts(10, fluxXsecPenalties);
    systlist.insert(systlist.end(), fluxlist.begin(), fluxlist.end());
  }
  if (detsyst) {
    if (useND) systlist.insert(systlist.end(), detlist_nd.begin(), detlist_nd.end());
    switch(det){
    case kFDApp: systlist.insert(systlist.end(), detlist_app.begin(), detlist_app.end()); break;
    case kFDDis: systlist.insert(systlist.end(), detlist_dis.begin(), detlist_dis.end()); break;
    case kFDAll: systlist.insert(systlist.end(), detlist_fd.begin(), detlist_fd.end()); break;
    }
  }
  if (xsecsyst) {
    std::vector<const ISyst*> xseclist = GetGenieSysts(GetGenieWeightNames(), fluxXsecPenalties);    
    systlist.insert(systlist.end(), xseclist.begin(), xseclist.end());
  }

  return systlist;
};

void RemoveSysts(std::vector<const ISyst *> &systlist,
                 std::vector<std::string> const &namesToRemove) {
  systlist.erase(std::remove_if(systlist.begin(), systlist.end(),
                                [&](const ISyst *s) {
                                  return (std::find(namesToRemove.begin(),
                                                    namesToRemove.end(),
                                                    s->ShortName()) !=
                                          namesToRemove.end());
                                }),
                 systlist.end());
}

TH2D* make_corr_from_covar(TH2D* covar){

  TH2D *corr = (TH2D*)covar->Clone();
  corr      ->SetNameTitle("corr", "corr");

  for (int i = 0; i < covar->GetNbinsX(); ++i){
    double istddev = sqrt(covar->GetBinContent(i+1, i+1));
    for (int j = 0; j < covar->GetNbinsX(); ++j){
      double jstddev  = sqrt(covar->GetBinContent(j+1, j+1));
      double new_corr = covar->GetBinContent(i+1, j+1)/istddev/jstddev; 
      corr ->SetBinContent(i+1, j+1, new_corr);
    }
  }
  return corr;
}

// Wow, this is ugly
std::vector<unique_ptr<PredictionInterp> > GetPredictionInterps(std::string fileName,
								std::vector<const ISyst*> systlist,
								bool reload=false){

  if(reload || TFile(fileName.c_str()).IsZombie()){
    
    // Need to start some loaders
    Loaders FD_FHC_loaders;
    Loaders FD_RHC_loaders;
    
    // Now fill the loaders
    FD_FHC_loaders .AddLoader(&FD_loaderFHCNumu,  caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNonSwap);
    FD_FHC_loaders .AddLoader(&FD_loaderFHCNue,   caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNueSwap);
    FD_FHC_loaders .AddLoader(&FD_loaderFHCNutau, caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNuTauSwap);

    FD_RHC_loaders .AddLoader(&FD_loaderRHCNumu,  caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNonSwap);
    FD_RHC_loaders .AddLoader(&FD_loaderRHCNue,   caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNueSwap);
    FD_RHC_loaders .AddLoader(&FD_loaderRHCNutau, caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNuTauSwap);

    NoExtrapPredictionGenerator genFDNumuFHC(axRecoEnuFDnumu, kPassFD_CVN_NUMU && kIsTrueFV, kGENIEWeights);

    NoExtrapPredictionGenerator genFDNumuRHC(axRecoEnuFDnumu, kPassFD_CVN_NUMU && kIsTrueFV, kGENIEWeights);
    
    NoExtrapPredictionGenerator genFDNueFHC(axRecoEnuFDnue, kPassFD_CVN_NUE && kIsTrueFV, kGENIEWeights);
    
    NoExtrapPredictionGenerator genFDNueRHC(axRecoEnuFDnue, kPassFD_CVN_NUE && kIsTrueFV, kGENIEWeights);
  
    // CW: Still need loaders at this stage for ND...
    NoOscPredictionGenerator genNDNumuFHC(ND_loaderFHC, axErecYrecND, kPassND_FHC_NUMU && kIsTrueFV, kGENIEWeights);
    
    NoOscPredictionGenerator genNDNumuRHC(ND_loaderRHC, axErecYrecND, kPassND_RHC_NUMU && kIsTrueFV, kGENIEWeights);
    
    osc::IOscCalculatorAdjustable* this_calc = NuFitOscCalc(1);      
    
    PredictionInterp predInterpFDNumuFHC(systlist,
					 this_calc, 
					 genFDNumuFHC, 
					 FD_FHC_loaders);
    
    PredictionInterp predInterpFDNumuRHC(systlist,
					 this_calc, 
					 genFDNumuRHC, 
					 FD_RHC_loaders);
    
    PredictionInterp predInterpFDNueFHC(systlist,
					this_calc, 
					genFDNueFHC, 
					FD_FHC_loaders);
      
    PredictionInterp predInterpFDNueRHC(systlist,
					this_calc, 
					genFDNueRHC, 
					FD_RHC_loaders);

    PredictionInterp predInterpNDNumuFHC(systlist,
					 this_calc,
					 genNDNumuFHC,
					 dummyLoaders);

    PredictionInterp predInterpNDNumuRHC(systlist,
					 this_calc,
					 genNDNumuRHC,
					 dummyLoaders);

    // Start all of the loaders
    ND_loaderFHC  .Go();
    ND_loaderRHC  .Go();
    FD_FHC_loaders.Go();
    FD_RHC_loaders.Go();
    
    TFile fout(fileName.c_str(), "RECREATE");
    std::cout << "Saving FD FHC numu" << std::endl;
    predInterpFDNumuFHC.SaveTo(fout.mkdir("fd_interp_numu_fhc"));
    std::cout << "Saving FD FHC nue" << std::endl;
    predInterpFDNueFHC .SaveTo(fout.mkdir("fd_interp_nue_fhc"));
    std::cout << "Saving FD RHC numu" << std::endl;
    predInterpFDNumuRHC.SaveTo(fout.mkdir("fd_interp_numu_rhc"));
    std::cout << "Saving FD RHC nue" << std::endl;
    predInterpFDNueRHC .SaveTo(fout.mkdir("fd_interp_nue_rhc"));
    std::cout << "Saving ND FHC" << std::endl;
    predInterpNDNumuFHC.SaveTo(fout.mkdir("nd_interp_numu_fhc"));
    std::cout << "Saving ND RHC" << std::endl;
    predInterpNDNumuRHC.SaveTo(fout.mkdir("nd_interp_numu_rhc"));
    fout .Close();
  }
  // Argh so ugly
  std::vector<unique_ptr<PredictionInterp> > return_list;
  TFile fin(fileName.c_str());
  std::cout << "Retrieving FD FHC numu" << std::endl;
  return_list.push_back(std::unique_ptr<ana::PredictionInterp> (LoadFrom<PredictionInterp>(fin.GetDirectory("fd_interp_numu_fhc")).release()));
  std::cout << "Retrieving FD FHC nue" << std::endl;
  return_list.push_back(std::unique_ptr<ana::PredictionInterp> (LoadFrom<PredictionInterp>(fin.GetDirectory("fd_interp_nue_fhc")).release()));
  std::cout << "Retrieving FD RHC numu" << std::endl;
  return_list.push_back(std::unique_ptr<ana::PredictionInterp> (LoadFrom<PredictionInterp>(fin.GetDirectory("fd_interp_numu_rhc")).release()));
  std::cout << "Retrieving FD RHC nue" << std::endl;
  return_list.push_back(std::unique_ptr<ana::PredictionInterp> (LoadFrom<PredictionInterp>(fin.GetDirectory("fd_interp_nue_rhc")).release()));
  std::cout << "Retrieving ND FHC" << std::endl;
  return_list.push_back(std::unique_ptr<ana::PredictionInterp> (LoadFrom<PredictionInterp>(fin.GetDirectory("nd_interp_numu_fhc")).release()));
  std::cout << "Retrieving ND RHC" << std::endl;
  return_list.push_back(std::unique_ptr<ana::PredictionInterp> (LoadFrom<PredictionInterp>(fin.GetDirectory("nd_interp_numu_rhc")).release()));
  fin.Close();
  return return_list;
};
