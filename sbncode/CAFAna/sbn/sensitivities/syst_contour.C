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

#include <vector>

const double sbndPOT = kPOTnominal;
const double icarusPOT = kPOTnominal;
const double uboonePOT = 1.3e21;

void syst_contour()
{
  //  GetSBNWeightSysts(); // initialize
  const std::vector<const ISyst*>& systs = GetSBNWeightSysts();

  TFile fin("cafe_state_syst_numu.root");

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

  osc::NoOscillations unosc;
  for(const ISyst* s: systs){
    p_nd->DebugPlot(s, &unosc);
    gPad->Print(TString::Format("plots/debug_nd_%s.pdf", s->ShortName().c_str()).Data());
    p_fd->DebugPlot(s, &unosc);
    gPad->Print(TString::Format("plots/debug_fd_%s.pdf", s->ShortName().c_str()).Data());


    p_nd->Predict(&unosc).ToTH1(sbndPOT)->Draw("hist");
    p_nd->PredictSyst(&unosc, SystShifts(s, -1)).ToTH1(sbndPOT, kBlue)->Draw("hist same");
    p_nd->PredictSyst(&unosc, SystShifts(s, +1)).ToTH1(sbndPOT, kRed)->Draw("hist same");
    p_nd->Predict(&unosc).ToTH1(sbndPOT)->Draw("hist same");
    leg_updn->Draw();
    gPad->Print(TString::Format("plots/spect_nd_%s.pdf", s->ShortName().c_str()).Data());

    p_fd->Predict(&unosc).ToTH1(icarusPOT)->Draw("hist");
    p_fd->PredictSyst(&unosc, SystShifts(s, -1)).ToTH1(icarusPOT, kBlue)->Draw("hist same");
    p_fd->PredictSyst(&unosc, SystShifts(s, +1)).ToTH1(icarusPOT, kRed)->Draw("hist same");
    p_fd->Predict(&unosc).ToTH1(icarusPOT)->Draw("hist same");
    leg_updn->Draw();
    gPad->Print(TString::Format("plots/spect_fd_%s.pdf", s->ShortName().c_str()).Data());

    p_ub->Predict(&unosc).ToTH1(uboonePOT)->Draw("hist");
    p_ub->PredictSyst(&unosc, SystShifts(s, -1)).ToTH1(icarusPOT, kBlue)->Draw("hist same");
    p_ub->PredictSyst(&unosc, SystShifts(s, +1)).ToTH1(icarusPOT, kRed)->Draw("hist same");
    p_ub->Predict(&unosc).ToTH1(icarusPOT)->Draw("hist same");
    leg_updn->Draw();
    gPad->Print(TString::Format("plots/spect_ub_%s.pdf", s->ShortName().c_str()).Data());

    if(fabs(p_fd->PredictSyst(&unosc, SystShifts(s, +1)).Integral(1e20)/p_fd->Predict(&unosc).Integral(1e20)-1) > .01) bigsysts.push_back(s);
  }

  std::cout << bigsysts.size() << " big systs out of " << systs.size() << std::endl;
  for(const ISyst* s: bigsysts) std::cout << s->ShortName() << " ";
  std::cout << std::endl;


  TLegend* leg = new TLegend(.15, .15, .45, .4);
  leg->SetFillStyle(0);
  dummy->SetLineColor(kRed);
  leg->AddEntry(dummy->Clone(), "SBND only", "l");
  dummy->SetLineColor(kGreen+2);
  leg->AddEntry(dummy->Clone(), "Microboone only", "l");
  dummy->SetLineColor(kBlue);
  leg->AddEntry(dummy->Clone(), "Icarus only", "l");
  dummy->SetLineColor(kBlack);
  leg->AddEntry(dummy->Clone(), "Combined fit", "l");
  dummy->SetLineStyle(7);
  leg->AddEntry(dummy->Clone(), "Stats only", "l");
  leg->Draw();


  OscCalcSterileApproxAdjustable* calc = DefaultSterileApproxCalc();

  //Define fit axes
  const FitAxis kAxSinSq2ThetaMuMu(&kFitSinSq2ThetaMuMu, 20/*40*/, 1e-3, 1, true);
  const FitAxis kAxDmSq(&kFitDmSqSterile, 40, 1e-2, 1e2, true);

  // We'll call zero nominal
  calc->SetL(kBaselineSBND);
  const Spectrum data_nd = p_nd->Predict(calc).FakeData(sbndPOT);
  calc->SetL(kBaselineIcarus);
  const Spectrum data_fd = p_fd->Predict(calc).FakeData(icarusPOT);
  calc->SetL(kBaselineMicroBoone);
  const Spectrum data_ub = p_ub->Predict(calc).FakeData(uboonePOT);

  SingleSampleExperiment expt_nd(p_nd, data_nd);
  SingleSampleExperiment expt_fd(p_fd, data_fd);
  SingleSampleExperiment expt_ub(p_ub, data_ub);

  MultiExperimentSBN multiExpt({&expt_nd, &expt_fd, &expt_ub}, {kSBND, kICARUS, kMicroBoone});


  Surface surf_nom(&multiExpt, calc, kAxSinSq2ThetaMuMu, kAxDmSq);
  calc->SetL(kBaselineSBND);
  Surface surf_nom_nd(&expt_nd, calc, kAxSinSq2ThetaMuMu, kAxDmSq);
  calc->SetL(kBaselineIcarus);
  Surface surf_nom_fd(&expt_fd, calc, kAxSinSq2ThetaMuMu, kAxDmSq);
  calc->SetL(kBaselineMicroBoone);
  Surface surf_nom_ub(&expt_ub, calc, kAxSinSq2ThetaMuMu, kAxDmSq);
    

  std::vector<std::vector<const ISyst*>> slists;
  slists.push_back(bigsysts);
  for(const ISyst* s: systs) slists.emplace_back(1, s); // and then each

  for(const std::vector<const ISyst*> slist: slists){
    new TCanvas;

    Surface surf_syst(&multiExpt, calc, kAxSinSq2ThetaMuMu, kAxDmSq, {}, slist);

    calc->SetL(kBaselineSBND);
    Surface surf_syst_nd(&expt_nd, calc, kAxSinSq2ThetaMuMu, kAxDmSq, {}, slist);
    calc->SetL(kBaselineIcarus);
    Surface surf_syst_fd(&expt_fd, calc, kAxSinSq2ThetaMuMu, kAxDmSq, {}, slist);
    calc->SetL(kBaselineMicroBoone);
    Surface surf_syst_ub(&expt_ub, calc, kAxSinSq2ThetaMuMu, kAxDmSq, {}, slist);


    //TH2* crit90 = Gaussian90Percent1D1Sided(surf_nom);
    TH2* crit3sig = Gaussian3Sigma1D1Sided(surf_nom);
    //TH2* crit5sig = Gaussian5Sigma1D1Sided(surf_nom);

    std::string suffix = "big";
    if(slist.size() == 1) suffix = slist[0]->ShortName();

    surf_nom_nd.SetTitle((suffix + " - 3#sigma C.L.").c_str());

    surf_nom_nd.DrawContour(crit3sig, 7, kRed);
    surf_nom_fd.DrawContour(crit3sig, 7, kBlue);
    surf_nom_ub.DrawContour(crit3sig, 7, kGreen+2);
    surf_nom.DrawContour(crit3sig, 7, kBlack);

    surf_syst_nd.DrawContour(crit3sig, kSolid, kRed);
    surf_syst_fd.DrawContour(crit3sig, kSolid, kBlue);
    surf_syst_ub.DrawContour(crit3sig, kSolid, kGreen+2);
    surf_syst.DrawContour(crit3sig, kSolid, kBlack);

    leg->Draw();
 
    gPad->Print(TString::Format("plots/conts_%s.pdf", suffix.c_str()).Data());

    if(slist.size() == 1){
      new TCanvas;
      TH2* h = surf_syst.GetProfiledHists()[0];
      h->SetMinimum(-2);
      h->SetMaximum(+2);
      h->Draw("colz");
      gPad->SetLogx();
      gPad->SetLogy();
      gPad->Print(TString::Format("plots/prof_%s.pdf", suffix.c_str()).Data());
    }
  } // end for s
}
