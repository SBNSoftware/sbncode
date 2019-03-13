#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Analysis/Calcs.h"
#include "CAFAna/Analysis/Plots.h"
#include "CAFAna/Analysis/Style.h"
#include "CAFAna/Analysis/Surface.h"
#include "CAFAna/Analysis/Fit.h"
#include "CAFAna/Vars/FitVars.h"
#include "CAFAna/Experiment/SingleSampleExperiment.h"
#include "CAFAna/Experiment/MultiExperiment.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Systs/DUNEXSecSysts.h"

using namespace ana;

#include "StandardRecord/StandardRecord.h"
#include "OscLib/func/IOscCalculator.h"

#include "Utilities/rootlogon.C"

#include "TCanvas.h"
#include "TGraph.h"
#include "TH2.h"
#include "TLegend.h"
#include "TPad.h"

void Legend()
{
  TLegend* leg = new TLegend(.6, .6, .9, .85);
  leg->SetFillStyle(0);

  TH1* dummy = new TH1F("", "", 1, 0, 1);
  dummy->SetMarkerStyle(kFullCircle);
  leg->AddEntry(dummy->Clone(), "Fake Data", "lep");
  dummy->SetLineColor(kTotalMCColor);
  leg->AddEntry(dummy->Clone(), "Total MC", "l");
  dummy->SetLineColor(kNCBackgroundColor);
  leg->AddEntry(dummy->Clone(), "NC", "l");
  dummy->SetLineColor(kNumuBackgroundColor);
  leg->AddEntry(dummy->Clone(), "#nu_{#mu} CC", "l");
  dummy->SetLineColor(kBeamNueBackgroundColor);
  leg->AddEntry(dummy->Clone(), "Beam #nu_{e} CC", "l");

  leg->Draw();
}

void test_dune()
{
  rootlogon(); // style

  // POT/yr * 3.5yrs * mass correction
  const double pot = 3.5 * 1.47e21 * 40/1.13;

  SpectrumLoader loaderNonswap("/dune/data/users/marshalc/CAFs/mcc10_test/FarDetector/FHC_nonswap.root");
  SpectrumLoader loaderNue("/dune/data/users/marshalc/CAFs/mcc10_test/FarDetector/FHC_nueswap.root");
  SpectrumLoaderBase& loaderNuTau = kNullLoader;

  const Var Enu_recoNumu = SIMPLEVAR(dune.Ev_reco_numu);
  const Var Enu_recoNue = SIMPLEVAR(dune.Ev_reco_nue);
  const Cut kSelNumu = SIMPLEVAR(dune.cvnnumu) > 0.7;
  const Cut kSelNue = SIMPLEVAR(dune.cvnnue) > 0.7;

  PredictionNoExtrap predNumuPID(loaderNonswap, loaderNue, loaderNuTau, "PID", Binning::Simple(60, -.1, +1.1), SIMPLEVAR(dune.cvnnumu), kNoCut);

  PredictionNoExtrap predNuePID(loaderNonswap, loaderNue, loaderNuTau, "PID", Binning::Simple(60, -.1, +1.1), SIMPLEVAR(dune.cvnnue), kNoCut);

  PredictionNoExtrap pred(loaderNonswap, loaderNue, loaderNuTau, "Reconstructed E (GeV)", Binning::Simple(80, 0, 10), Enu_recoNumu, kSelNumu);

  PredictionNoExtrap predNue(loaderNonswap, loaderNue, loaderNuTau, "Reconstructed E (GeV)", Binning::Simple(24, 0, 6), Enu_recoNue, kSelNue);

  /*
  // separate by true interaction category
  std::vector<Cut> truthcuts;
  for( int i = 0; i < 32; ++i ) {
    truthcuts.push_back( kVALORCategory == i );
  }
  const HistAxis axis("Reconstructed energy (GeV)",
                      Binning::Simple(40, 0, 10),
                      Enu_recoNumu);
  PredictionScaleComp predNumu(loaderNonswap, loaderNue, loaderNuTau, axis, kSelNumu, truthcuts);
  */

  loaderNonswap.Go();
  loaderNue.Go();

  SaveToFile(pred, "pred_numu.root", "pred");
  SaveToFile(predNue, "pred_nue.root", "pred");

  osc::IOscCalculatorAdjustable* calc = DefaultOscCalc();
  calc->SetL(1300);
  calc->SetdCP(TMath::Pi()*1.5);


  // Standard DUNE numbers from Matt Bass
  calc->SetTh12(0.5857);
  calc->SetTh13(0.148);
  calc->SetTh23(0.726);
  calc->SetDmsq21(0.000075);
  calc->SetDmsq32(0.002524-0.000075); // quoted value is 31

  // One sigma errors
  // (t12,t13,t23,dm21,dm32)=(0.023,0.018,0.058,0.0,0.024,0.016)


  Spectrum mock = pred.Predict(calc).FakeData(pot);
  SingleSampleExperiment expt(&pred, mock);

  Spectrum mockNue = predNue.Predict(calc).FakeData(pot);
  SingleSampleExperiment exptNue(&predNue, mockNue);

  Spectrum mockNuePID = predNuePID.Predict(calc).FakeData(pot);
  Spectrum mockNumuPID = predNumuPID.Predict(calc).FakeData(pot);

  Surface surf(&expt, calc, &kFitSinSqTheta23, 20, .4, .6, &kFitDmSq32Scaled, 20, 2.35, 2.5);

  MultiExperiment me({&expt, &exptNue});

  Surface surfNue(&me, calc, &kFitDeltaInPiUnits, 20, 0, 2, &kFitSinSqTheta23, 20, .4, .6);


  new TCanvas;
  TH1* h3 = DataMCComparisonComponents(mock, &pred, calc);
  h3->SetTitle("#nu_{#mu} selection (CVNm>0.7) 3.5yrs #times 40kt");
  CenterTitles(h3);
  Legend();
  gPad->Print("components.eps");
  gPad->Print("components.C");

  new TCanvas;
  TH1* h4 = DataMCComparisonComponents(mockNue, &predNue, calc);
  h4->SetTitle("#nu_{e} selection (CVNe>0.7) 3.5yrs #times 40kt");
  CenterTitles(h4);
  Legend();
  gPad->Print("components_nue.eps");
  gPad->Print("components_nue.C");

  new TCanvas;
  TH1* h2 = DataMCComparisonComponents(mockNumuPID, &predNumuPID, calc);
  h2->SetTitle("#nu_{#mu} selection (CVNm>0.7) 3.5yrs #times 40kt");
  CenterTitles(h2);
  Legend();
  h2->GetYaxis()->SetRangeUser(1, 1e4);
  gPad->SetLogy();
  gPad->Print("components_pid.eps");
  gPad->Print("components_pid.C");

  new TCanvas;
  TH1* h = DataMCComparisonComponents(mockNuePID, &predNuePID, calc);
  h->SetTitle("#nu_{e} selection (CVNe>0.7) 3.5yrs #times 40kt");
  CenterTitles(h);
  Legend();
  h->GetYaxis()->SetRangeUser(0, 600);
  gPad->Print("components_nue_pid.eps");
  gPad->Print("components_nue_pid.C");

  new TCanvas;
  surf.DrawContour(Gaussian68Percent2D(surf), 7, kRed);
  surf.DrawContour(Gaussian90Percent2D(surf), kSolid, kRed);
  surf.DrawBestFit(kRed);
  gPad->Print("cont.eps");

  new TCanvas;
  surfNue.DrawContour(Gaussian68Percent2D(surfNue), 7, kRed);
  surfNue.DrawContour(Gaussian90Percent2D(surfNue), kSolid, kRed);
  surfNue.DrawBestFit(kRed);
  gPad->Print("cont_nue.eps");


  new TCanvas;

  // This is a very cheesy way to make the McD plot - would have to be very
  // different if we were varying any other parameters
  calc->SetdCP(0);
  Spectrum zeroNumu = pred.Predict(calc).FakeData(pot);
  Spectrum zeroNue = predNue.Predict(calc).FakeData(pot);
  calc->SetdCP(TMath::Pi());
  Spectrum oneNumu = pred.Predict(calc).FakeData(pot);
  Spectrum oneNue = predNue.Predict(calc).FakeData(pot);

  TGraph* g = new TGraph;

  for(int i = -100; i <= 200; ++i){
    calc->SetdCP(i/100.*TMath::Pi());

    Spectrum mockNumu = pred.Predict(calc).FakeData(pot);
    Spectrum mockNue = predNue.Predict(calc).FakeData(pot);

    const double llZero = //LogLikelihood(zeroNumu.ToTH1(pot), mockNumu.ToTH1(pot))+
      LogLikelihood(zeroNue.ToTH1(pot), mockNue.ToTH1(pot));

    const double llOne = //LogLikelihood(oneNumu.ToTH1(pot), mockNumu.ToTH1(pot))+
      LogLikelihood(oneNue.ToTH1(pot), mockNue.ToTH1(pot));

    const double ll = std::min(llZero, llOne);

    g->SetPoint(g->GetN(), i/100., sqrt(ll));
  }

  TH2* axes = new TH2F("", "3.5yrs #times 40kt, stats only;#delta / #pi;#sqrt{#chi^{2}}", 100, -1, +1, 100, 0, 6);
  CenterTitles(axes);
  axes->Draw();
  g->SetLineWidth(2);
  g->Draw("l same");

  gPad->SetGridx();
  gPad->SetGridy();

  gPad->Print("mcd.eps");
  gPad->Print("mcd.C");


  /*
  SystShifts fakeDataShift(predNumu.GetSysts()[k_int_nu_MEC_dummy], +1);
  Spectrum fake = predNumu.PredictSyst(calc, fakeDataShift).FakeData(1.47e21);

  fake.ToTH1(1.47e21, kRed)->Draw("hist");
  predNumu.Predict(calc).ToTH1(1.47e21)->Draw("hist same");

  SingleSampleExperiment mecexpt(&predNumu, fake);

  Fitter fit(&mecexpt, {}, predNumu.GetSysts());
  SystShifts seed = SystShifts::Nominal();
  fit.Fit(seed);

  SystShifts bfs(predNumu.GetSysts()[k_int_nu_MEC_dummy], seed.GetShift(predNumu.GetSysts()[k_int_nu_MEC_dummy]));
  Spectrum bf = predNumu.PredictSyst(calc, bfs);
  Spectrum bf2 = predNumu.PredictSyst(calc, seed);
  bf.ToTH1(1.47e21, kBlue)->Draw("hist same");
  bf2.ToTH1(1.47e21, kBlue, 7)->Draw("hist same");

  gPad->Print("FD_mec_fit.eps");
  */
}
