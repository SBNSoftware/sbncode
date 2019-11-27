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

void SpectrumPlots(const PredictionInterp* p, const ISyst* s, double pot, const char* detName)
{
  osc::NoOscillations unosc;

  p->DebugPlot(s, &unosc);
  gPad->Print(TString::Format("plots/debug_%s_%s.pdf", detName, s->ShortName().c_str()).Data());

  p->Predict(&unosc).ToTH1(pot)->Draw("hist");
  p->PredictSyst(&unosc, SystShifts(s, -1)).ToTH1(pot, kBlue)->Draw("hist same");
  p->PredictSyst(&unosc, SystShifts(s, +1)).ToTH1(pot, kRed )->Draw("hist same");
  p->Predict(&unosc).ToTH1(pot)->Draw("hist same");

  TLegend* leg = new TLegend(.6, .6, .85, .85);
  leg->SetFillStyle(0);
  TH1* dummy = new TH1F("", "", 1, 0, 1);
  leg->AddEntry(dummy->Clone(), "Nominal", "l");
  dummy->SetLineColor(kBlue);
  leg->AddEntry(dummy->Clone(), "-1#sigma", "l");
  dummy->SetLineColor(kRed);
  leg->AddEntry(dummy->Clone(), "+1#sigma", "l");

  leg->Draw();

  gPad->Print(TString::Format("plots/spect_%s_%s.pdf", detName, s->ShortName().c_str()).Data());
}


void syst_contour_nue()
{
  const std::vector<const ISyst*>& systs = GetSBNWeightSysts();

  TFile fin("cafe_state_syst_nue.root");

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

  // We'll accumulate all the systematics that have more than a 1% effect on
  // the normalization in this vector.
  std::vector<const ISyst*> bigsysts;

  osc::NoOscillations unosc;
  for(const ISyst* s: systs){
    SpectrumPlots(p_nd, s, sbndPOT, "nd");
    SpectrumPlots(p_fd, s, icarusPOT, "fd");
    SpectrumPlots(p_ub, s, uboonePOT, "ub");

    osc::NoOscillations unosc;
    if(fabs(p_fd->PredictSyst(&unosc, SystShifts(s, +1)).Integral(1e20)/p_fd->Predict(&unosc).Integral(1e20)-1) > .01) bigsysts.push_back(s);
  } // end for s

  std::cout << bigsysts.size() << " big systs out of " << systs.size() << std::endl;
  for(const ISyst* s: bigsysts) std::cout << s->ShortName() << " ";
  std::cout << std::endl;


  OscCalcSterileApproxAdjustable* calc = DefaultSterileApproxCalc();

  //Define fit axes
  const FitAxis kAxSinSq2ThetaMuE(&kFitSinSq2ThetaMuE, 20/*40*/, 1e-4, 1, true);
  const FitAxis kAxDmSq(&kFitDmSqSterile, 40, 1e-2, 1e2, true);

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


  Surface surf_nom(&multiExpt, calc, kAxSinSq2ThetaMuE, kAxDmSq);
  calc->SetL(kBaselineSBND);
  Surface surf_nom_nd(&expt_nd, calc, kAxSinSq2ThetaMuE, kAxDmSq);
  calc->SetL(kBaselineIcarus);
  Surface surf_nom_fd(&expt_fd, calc, kAxSinSq2ThetaMuE, kAxDmSq);
  calc->SetL(kBaselineMicroBoone);
  Surface surf_nom_ub(&expt_ub, calc, kAxSinSq2ThetaMuE, kAxDmSq);


  //TH2* crit90 = Gaussian90Percent1D1Sided(surf_nom);
  TH2* crit3sig = Gaussian3Sigma1D1Sided(surf_nom);
  //TH2* crit5sig = Gaussian5Sigma1D1Sided(surf_nom);

  new TCanvas;
  surf_nom_nd.DrawContour(crit3sig, kSolid, kRed);
  surf_nom_fd.DrawContour(crit3sig, kSolid, kBlue);
  surf_nom_ub.DrawContour(crit3sig, kSolid, kGreen+2);
  surf_nom.DrawContour(crit3sig, kSolid, kBlack);


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


  leg->Draw();

  gPad->Print("plots/conts_nom.pdf");

  dummy->SetLineStyle(7);
  leg->AddEntry(dummy->Clone(), "Stats only", "l");

  // We'll make contours first for all the big systematics together, and then
  // we'll make the effect of each systematic independently.
  std::vector<std::vector<const ISyst*>> slists;
  slists.push_back(bigsysts); // all systs
  for(const ISyst* s: systs) slists.emplace_back(1, s); // then each individually

  for(const std::vector<const ISyst*> slist: slists){
    new TCanvas;

    Surface surf_syst(&multiExpt, calc, kAxSinSq2ThetaMuE, kAxDmSq, {}, slist);

    calc->SetL(kBaselineSBND);
    Surface surf_syst_nd(&expt_nd, calc, kAxSinSq2ThetaMuE, kAxDmSq, {}, slist);
    calc->SetL(kBaselineIcarus);
    Surface surf_syst_fd(&expt_fd, calc, kAxSinSq2ThetaMuE, kAxDmSq, {}, slist);
    calc->SetL(kBaselineMicroBoone);
    Surface surf_syst_ub(&expt_ub, calc, kAxSinSq2ThetaMuE, kAxDmSq, {}, slist);

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

    // Show the profiled value of each syst for the fits that use only a single
    // syst
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
