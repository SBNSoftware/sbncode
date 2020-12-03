#include "CAFAna/Core/LoadFromFile.h"
#include "CAFAna/Core/OscCalcSterileApprox.h"
#include "CAFAna/Vars/FitVarsSterileApprox.h"
#include "CAFAna/Prediction/PredictionInterp.h"
#include "CAFAna/Experiment/SingleSampleExperiment.h"
#include "CAFAna/Experiment/MultiExperimentSBN.h"
#include "CAFAna/Experiment/CountingExperiment.h"
#include "CAFAna/Analysis/ExpInfo.h"
#include "CAFAna/Analysis/Surface.h"
#include "CAFAna/Systs/SBNWeightSysts.h"
#include "CAFAna/Systs/SBNFluxSysts.h"
#include "CAFAna/Systs/Systs.h"
using namespace ana;

#include "OscLib/IOscCalculator.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH2.h"
#include "TLegend.h"
#include "TString.h"
#include "TRandom.h"

#include <vector>

const double sbndPOT = kPOTnominal;
const double icarusPOT = kPOTnominal;
const double uboonePOT = 1.3e21;

const std::string numuStr = "numu";
const std::string nueStr = "nue";

void syst_spec(const std::string anatype = numuStr)
{
  const std::vector<const ISyst*>& systs = GetSBNWeightSysts();

  auto systs_flux = GetSBNFluxWeightSysts();
  auto systs_genie = GetSBNGenieWeightSysts();
  
  GetSBNFluxHadronSysts(30);
  GetMECSyst();

  std::vector<const ISyst*> systs_to_process;

  std::vector<std::string> syst_names{"expskin_FluxUnisim","horncurrent_FluxUnisim","kminus_PrimaryHadronNormalization","kplus_PrimaryHadronFeynmanScaling","kzero_PrimaryHadronSanfordWang","nucleoninexsec_FluxUnisim","nucleonqexsec_FluxUnisim","nucleontotxsec_FluxUnisim","piminus_PrimaryHadronSWCentralSplineVariation","pioninexsec_FluxUnisim","pionqexsec_FluxUnisim","piontotxsec_FluxUnisim","piplus_PrimaryHadronSWCentralSplineVariation","genie_ccresAxial_Genie","genie_ncresAxial_Genie","genie_qema_Genie","genie_NC_Genie","genie_NonResRvbarp1pi_Genie","genie_NonResRvbarp2pi_Genie","genie_NonResRvp1pi_Genie","genie_NonResRvp2pi_Genie","genie_NonResRvbarp1piAlt_Genie","genie_NonResRvbarp2piAlt_Genie","genie_NonResRvp1piAlt_Genie","genie_NonResRvp2piAlt_Genie"};

  for (auto s : systs) {
    for (auto n : syst_names) if (n == s->ShortName()) systs_to_process.push_back(s);
  }

  for(const ISyst* s: GetSBNFluxHadronSysts(30)) systs_to_process.push_back(s);
  
  std::map<std::string, std::vector<const ISyst*>> slists;
  slists["prop_systs"] = systs_to_process;
  slists["flux_noPCA"] = systs_flux;
  slists["genie"] = systs_genie;

  const char* name_in;
  const char* name_out;
  if (anatype == numuStr) {
    name_in = "cafe_state_syst_numu.root";
    name_out = "prefit_syst_spec_numu.root";
  }
  else if (anatype == nueStr) {
    name_in = "cafe_state_syst_nue.root";
    name_out = "prefit_syst_spec_nue.root";
  }
  else {
    std::cout << "Must specifiy nue or numu" << std::endl;
    return;
  }

  TFile fin(name_in);
  TFile fout(name_out, "RECREATE");

  PredictionInterp* p_nd = LoadFrom<PredictionInterp>(fin.GetDirectory("pred_nd")).release();
  PredictionInterp* p_fd = LoadFrom<PredictionInterp>(fin.GetDirectory("pred_fd")).release();
  PredictionInterp* p_ub = LoadFrom<PredictionInterp>(fin.GetDirectory("pred_ub")).release();
 
  osc::NoOscillations unosc;

  auto nom_nd = p_nd->Predict(&unosc).ToTH1(sbndPOT);
  auto nom_ub = p_ub->Predict(&unosc).ToTH1(uboonePOT);
  auto nom_fd = p_fd->Predict(&unosc).ToTH1(icarusPOT);
  nom_nd->Scale(1, "width");
  nom_ub->Scale(1, "width");
  nom_fd->Scale(1, "width");
  nom_nd->Write("spec_nd_nominal");
  nom_ub->Write("spec_ub_nominal");
  nom_fd->Write("spec_fd_nominal");

  for(auto syst_pair: slists) {
    std::string suffix = syst_pair.first;
    std::vector<const ISyst*> slist = syst_pair.second;    

    std::vector<TH1*> hists;
    for (int j = 0; j < 1000; ++j){
      SystShifts shifts = SystShifts::RandomThrow(slist);
      hists.push_back(p_nd->PredictSyst(&unosc, shifts).ToTH1(sbndPOT));
    }
    double xbins[] = {0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.25, 1.5, 2., 2.5, 3.};
    TH1D *upper = new TH1D(("h1"+suffix).c_str(), "hist", 19, xbins);
    TH1D *lower = new TH1D(("h2"+suffix).c_str(), "hist", 19, xbins);
    for (int k = 1; k <= 19; ++k ) {
      std::vector<double> bincont;
      for (auto h : hists) bincont.push_back(h->GetBinContent(k));
      std::sort(bincont.begin(),bincont.end());
      upper->SetBinContent(k, bincont[840]);
      lower->SetBinContent(k, bincont[160]);
    }    
    upper->Scale(1, "width");
    lower->Scale(1, "width");
    upper->Write(("spec_nd_"+suffix+"_+1").c_str());
    lower->Write(("spec_nd_"+suffix+"_-1").c_str());

    std::vector<TH1*> hists2;
    for (int j = 0; j < 1000; ++j){
      SystShifts shifts = SystShifts::RandomThrow(slist);
      hists2.push_back(p_ub->PredictSyst(&unosc, shifts).ToTH1(uboonePOT));
    }
    TH1D *upper2 = new TH1D(("h3"+suffix).c_str(), "hist", 19, xbins);
    TH1D *lower2 = new TH1D(("h4"+suffix).c_str(), "hist", 19, xbins);
    for (int k = 1; k <= 19; ++k ) {
      std::vector<double> bincont;
      for (auto h : hists2) bincont.push_back(h->GetBinContent(k));
      std::sort(bincont.begin(),bincont.end());
      upper2->SetBinContent(k, bincont[840]);
      lower2->SetBinContent(k, bincont[160]);
    }    
    upper2->Scale(1, "width");
    lower2->Scale(1, "width");
    upper2->Write(("spec_ub_"+suffix+"_+1").c_str());
    lower2->Write(("spec_ub_"+suffix+"_-1").c_str());

    std::vector<TH1*> hists3;
    for (int j = 0; j < 1000; ++j){
      SystShifts shifts = SystShifts::RandomThrow(slist);
      hists3.push_back(p_fd->PredictSyst(&unosc, shifts).ToTH1(icarusPOT));
    }
    TH1D *upper3 = new TH1D(("h5"+suffix).c_str(), "hist", 19, xbins);
    TH1D *lower3 = new TH1D(("h6"+suffix).c_str(), "hist", 19, xbins);
    for (int k = 1; k <= 19; ++k ) {
      std::vector<double> bincont;
      for (auto h : hists3) bincont.push_back(h->GetBinContent(k));
      std::sort(bincont.begin(),bincont.end());
      upper3->SetBinContent(k, bincont[840]);
      lower3->SetBinContent(k, bincont[160]);
    }    
    upper3->Scale(1, "width");
    lower3->Scale(1, "width");
    upper3->Write(("spec_fd_"+suffix+"_+1").c_str());
    lower3->Write(("spec_fd_"+suffix+"_-1").c_str());
  }
}
