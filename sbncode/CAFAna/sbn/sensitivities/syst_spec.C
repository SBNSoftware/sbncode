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

void syst_spec()
{
  //  GetSBNWeightSysts(); // initialize
  const std::vector<const ISyst*>& systs = GetSBNWeightSysts();

  auto systs_flux = GetSBNFluxWeightSysts();
  auto systs_genie = GetSBNGenieWeightSysts();
  
  std::vector<std::vector<const ISyst*>> all_systs_vec;
  all_systs_vec.push_back(systs);
  all_systs_vec.push_back(systs_flux);
  all_systs_vec.push_back(systs_genie);
  
 std::string n[] = {"all", "flux", "genie"};

  TFile fin("cafe_state_syst_numu.root");
  TFile fout("syst_spec.root","RECREATE");


  PredictionInterp* p_nd = LoadFrom<PredictionInterp>(fin.GetDirectory("pred_nd")).release();
  PredictionInterp* p_fd = LoadFrom<PredictionInterp>(fin.GetDirectory("pred_fd")).release();
  PredictionInterp* p_ub = LoadFrom<PredictionInterp>(fin.GetDirectory("pred_ub")).release();

  TLegend* leg_updn = new TLegend(.6, .6, .85, .85);
  leg_updn->SetFillStyle(0);
  TH1* dummy = new TH1F("", "", 1, 0, 1);
  leg_updn->AddEntry(dummy->Clone(), "Nominal", "l");
  dummy->SetLineColor(kBlue);
  leg_updn->AddEntry(dummy->Clone(), "-1#sigma", "l");
  dummy->SetLineColor(kRed);
  leg_updn->AddEntry(dummy->Clone(), "+1#sigma", "l");

  std::vector<const ISyst*> bigsysts;

  // Some code from Chris' script syst_contour.C
  // I wasn't using it, but it can be used to look at individual systematics
  // So I left it, but it's commented out
  osc::NoOscillations unosc;
//  for(const ISyst* s: systs){
//    p_nd->DebugPlot(s, &unosc);
//    gPad->Print(TString::Format("plots/debug_nd_%s.pdf", s->ShortName().c_str()).Data());
//    p_fd->DebugPlot(s, &unosc);
//    gPad->Print(TString::Format("plots/debug_fd_%s.pdf", s->ShortName().c_str()).Data());
//
//
//    auto h3 = p_nd->PredictSyst(&unosc, SystShifts(s, +1)).ToTH1(sbndPOT, kRed);
//    h3->Write(("spect_nd_"+s->ShortName()+"_+1").c_str(), TObject::kOverwrite);
//    h3->Scale(1, "width");
//    h3->Draw("hist");
//    auto h1 = p_nd->Predict(&unosc).ToTH1(sbndPOT);
//    h1->Write(("spect_nd_"+s->ShortName()+"_no_shift").c_str(), TObject::kOverwrite);
//    h1->Scale(1, "width");
//    h1->Draw("hist same");
//    auto h2 = p_nd->PredictSyst(&unosc, SystShifts(s, -1)).ToTH1(sbndPOT, kBlue);
//    h2->Write(("spect_nd_"+s->ShortName()+"_-1").c_str(), TObject::kOverwrite);
//    h2->Scale(1, "width");
//    h2->Draw("hist same");
//    auto h4 = p_nd->Predict(&unosc).ToTH1(sbndPOT);
//    h4->Write(("spect_nd_"+s->ShortName()+"_no_shift2").c_str(), TObject::kOverwrite);
//    h4->Scale(1, "width");
//    h4->Draw("hist same");
//    leg_updn->Draw();
//    gPad->Print(TString::Format("spec_plots/spect_nd_%s.pdf", s->ShortName().c_str()).Data());
//
//    auto h5 = p_fd->Predict(&unosc).ToTH1(icarusPOT);
//    h5->Write(("spect_fd_"+s->ShortName()+"_no_shift").c_str(), TObject::kOverwrite);
//    h5->Scale(1, "width");
//    h5->Draw("hist");
//    auto h6 = p_fd->PredictSyst(&unosc, SystShifts(s, -1)).ToTH1(icarusPOT, kBlue);
//    h6->Write(("spect_fd_"+s->ShortName()+"_-1").c_str(), TObject::kOverwrite);
//    h6->Scale(1, "width");
//    h6->Draw("hist same");
//    auto h7 = p_fd->PredictSyst(&unosc, SystShifts(s, +1)).ToTH1(icarusPOT, kRed);
//    h7->Write(("spect_fd_"+s->ShortName()+"_+1").c_str(), TObject::kOverwrite);
//    h7->Scale(1, "width");
//    h7->Draw("hist same");
//    auto h8 = p_fd->Predict(&unosc).ToTH1(icarusPOT);
//    h8->Write(("spect_fd_"+s->ShortName()+"_no_shift2").c_str(), TObject::kOverwrite);
//    h8->Scale(1, "width");
//    h8->Draw("hist same");
//    leg_updn->Draw();
//    gPad->Print(TString::Format("spec_plots/spect_fd_%s.pdf", s->ShortName().c_str()).Data());
//
//    auto h9 = p_ub->Predict(&unosc).ToTH1(uboonePOT);
//    h9->Write(("spect_ub_"+s->ShortName()+"_no_shift").c_str(), TObject::kOverwrite);
//    h9->Scale(1, "width");
//    h9->Draw("hist");
//    auto h10 = p_ub->PredictSyst(&unosc, SystShifts(s, -1)).ToTH1(uboonePOT, kBlue);
//    h10->Write(("spect_ub_"+s->ShortName()+"_-1").c_str(), TObject::kOverwrite);
//    h10->Scale(1, "width");
//    h10->Draw("hist same");
//    auto h11 = p_ub->PredictSyst(&unosc, SystShifts(s, +1)).ToTH1(uboonePOT, kRed);
//    h11->Write(("spect_ub_"+s->ShortName()+"_+1").c_str(), TObject::kOverwrite);
//    h11->Scale(1, "width");
//    h11->Draw("hist same");
//    auto h12 = p_ub->Predict(&unosc).ToTH1(uboonePOT);
//    h12->Write(("spect_ub_"+s->ShortName()+"_no_shift2").c_str(), TObject::kOverwrite);
//    h12->Scale(1, "width");
//    h12->Draw("hist same");
//    leg_updn->Draw();
//    gPad->Print(TString::Format("spec_plots/spect_ub_%s.pdf", s->ShortName().c_str()).Data());
//
//    if(fabs(p_fd->PredictSyst(&unosc, SystShifts(s, +1)).Integral(1e20)/p_fd->Predict(&unosc).Integral(1e20)-1) > .01) bigsysts.push_back(s);
//  }
//
//  std::cout << bigsysts.size() << " big systs out of " << det_systs.size() << std::endl;
//  for(const ISyst* s: bigsysts) std::cout << s->ShortName() << " ";
//  std::cout << std::endl;


  for(int i = 0; i < 3; ++i) {
    SystShifts shifts;
    std::vector<TH1*> hists;
    for (int j = 0; j < 1000; ++j){
      for(auto s : all_systs_vec[i]) shifts.SetShift(s, gRandom->Gaus());
      hists.push_back(p_nd->PredictSyst(&unosc, shifts).ToTH1(sbndPOT));
    }
    double xbins[] = {0.2, 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1., 1.25, 1.5, 2., 2.5, 3.};
    TH1D *shifted = new TH1D(("h"+to_string(i)).c_str(), "hist",19,xbins);
    TH1D *lower = new TH1D("h2", "hist", 19, xbins);
    for (int k = 1; k <= 19; ++k ) {
      std::vector<double> bincont;
      for (auto h : hists) bincont.push_back(h->GetBinContent(k));
      std::sort(bincont.begin(),bincont.end());
      shifted->SetBinContent(k, bincont[840]);
      lower->SetBinContent(k, bincont[160]);
    }    
    
    shifted->Write(("spect_nd_"+n[i]+"_+1").c_str(), TObject::kOverwrite);
    lower->Write(("spect_nd_"+n[i]+"_-1").c_str(), TObject::kOverwrite);

    SystShifts shifts2;
    std::vector<TH1*> hists2;
    for (int j = 0; j < 1000; ++j){
      for(auto s : all_systs_vec[i]) shifts2.SetShift(s, gRandom->Gaus());
      hists2.push_back(p_fd->PredictSyst(&unosc, shifts2).ToTH1(icarusPOT));
    }
    TH1D *shifted2 = new TH1D(("h"+to_string(i)).c_str(), "hist",19,xbins);
    TH1D *lower2 = new TH1D("h2", "hist", 19, xbins);
    for (int k = 1; k <= 19; ++k ) {
      std::vector<double> bincont;
      for (auto h : hists2) bincont.push_back(h->GetBinContent(k));
      std::sort(bincont.begin(),bincont.end());
      shifted2->SetBinContent(k, bincont[840]);
      lower2->SetBinContent(k, bincont[160]);
    }    
    
    shifted2->Write(("spect_fd_"+n[i]+"_+1").c_str(), TObject::kOverwrite);
    lower2->Write(("spect_fd_"+n[i]+"_-1").c_str(), TObject::kOverwrite);

    SystShifts shifts3;
    std::vector<TH1*> hists3;
    for (int j = 0; j < 1000; ++j){
      for(auto s : all_systs_vec[i]) shifts3.SetShift(s, gRandom->Gaus());
      hists3.push_back(p_ub->PredictSyst(&unosc, shifts3).ToTH1(uboonePOT));
    }
    TH1D *shifted3 = new TH1D(("h"+to_string(i)).c_str(), "hist",19,xbins);
    TH1D *lower3 = new TH1D("h2", "hist", 19, xbins);
    for (int k = 1; k <= 19; ++k ) {
      std::vector<double> bincont;
      for (auto h : hists3) bincont.push_back(h->GetBinContent(k));
      std::sort(bincont.begin(),bincont.end());
      shifted3->SetBinContent(k, bincont[840]);
      lower3->SetBinContent(k, bincont[160]);
    }    
    
    shifted3->Write(("spect_ub_"+n[i]+"_+1").c_str(), TObject::kOverwrite);
    lower3->Write(("spect_ub_"+n[i]+"_-1").c_str(), TObject::kOverwrite);
  }
}
