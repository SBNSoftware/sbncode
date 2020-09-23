#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Prediction/PredictionSBNExtrap.h"
#include "CAFAna/Prediction/PredictionInterp.h"
#include "CAFAna/Analysis/Calcs.h"
#include "CAFAna/Analysis/FitAxis.h"

#include "CAFAna/Core/OscCalcSterileApprox.h"
#include "CAFAna/Vars/FitVarsSterileApprox.h"

#include "CAFAna/Analysis/Surface.h"
#include "CAFAna/Experiment/SingleSampleExperiment.h"
#include "CAFAna/Experiment/MultiExperimentSBN.h"
#include "CAFAna/Experiment/CountingExperiment.h"
#include "CAFAna/Analysis/ExpInfo.h"

#include "CAFAna/Systs/SBNWeightSysts.h"
#include "CAFAna/Systs/SBNFluxSysts.h"

#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"

using namespace ana;

// State file
const char* basicFname = "cafe_state_extrap_syst.root";

const double sbndPOT = kPOTnominal;
const double icarusPOT = kPOTnominal;

void nus_extrap(const char* stateFname = basicFname)
{
  if (TFile(stateFname).IsZombie()){
    std::cout << "Run make_state_extrap.C first!" << std::endl;
    return;
  }

  std::vector<const ISyst*> systs = GetSBNWeightSysts();
  for(const ISyst* s: GetSBNFluxHadronSysts(30)) systs.push_back(s);
//  systs.resize(10);

  std::cout << "Loading state from " << stateFname << std::endl; 
  TFile fin(stateFname);
  IPrediction* pred = ana::LoadFrom<IPrediction>(fin.GetDirectory("pred_syst")).release();
  IPrediction* predND = ana::LoadFrom<IPrediction>(fin.GetDirectory("predND_syst")).release();

  // Calculator
  OscCalcSterileApproxAdjustable* calc = DefaultSterileApproxCalc();
  // TODO - why is this necessary?
  ((PredictionInterp*)pred)->SetOscSeed(calc);

  // To make a fit we need to have a "data" spectrum to compare to our MC
  // Prediction object
  calc->SetL(kBaselineIcarus);
  const Spectrum data = pred->Predict(calc).FakeData(icarusPOT);
  SingleSampleExperiment expt(pred, data);

  calc->SetL(kBaselineSBND);
  const Spectrum dataND = predND->Predict(calc).FakeData(sbndPOT);
  CountingExperiment exptCount(predND, dataND);

  MultiExperimentSBN multiExpt({&expt, &exptCount}, {kICARUS, kSBND});

  //Define fit axes
  const FitAxis kAxSinSq2ThetaMuMu(&kFitSinSq2ThetaMuMu, 40, 1e-3, 1, true);
  const FitAxis kAxDmSq(&kFitDmSqSterile, 40, 2e-2, 1e2, true);

  const FitAxis kAxSinSq2ThetaMuMuCoarse(&kFitSinSq2ThetaMuMu, 8/*40*/, 1e-3, 1, true);
  const FitAxis kAxDmSqCoarse(&kFitDmSqSterile, 8/*40*/, 2e-2, 1e2, true);

  TCanvas* c1 = new TCanvas("c1");
  c1->SetLeftMargin(0.12);
  c1->SetBottomMargin(0.15);

  // A Surface evaluates the experiment's chisq across a grid
  Surface surf(&expt, calc,
               kAxSinSq2ThetaMuMuCoarse,
               kAxDmSqCoarse,
               {}, systs);

  Surface surfStat(&expt, calc,
                   kAxSinSq2ThetaMuMu,
                   kAxDmSq);

  calc->SetL(kBaselineSBND);
  Surface surfCount(&exptCount, calc,
                    kAxSinSq2ThetaMuMuCoarse,
                    kAxDmSqCoarse,
                    {}, systs);

  Surface surfCountStat(&exptCount, calc,
                        kAxSinSq2ThetaMuMu,
                        kAxDmSq);

  Surface surfJoint(&multiExpt, calc,
                    kAxSinSq2ThetaMuMuCoarse,
                    kAxDmSqCoarse,
                    {}, systs);

  Surface surfJointStat(&multiExpt, calc,
                        kAxSinSq2ThetaMuMu,
                        kAxDmSq);

  TH2* crit3sig = Gaussian3Sigma1D1Sided(surfStat);
  TH2* crit3sigCoarse = Gaussian3Sigma1D1Sided(surf);

  surfJointStat.DrawContour(crit3sig, 7, kRed);
  surfJoint.DrawContour(crit3sigCoarse, kSolid, kRed);

  surfStat.DrawContour(crit3sig, 7, kGreen+2);
  surf.DrawContour(crit3sigCoarse, kSolid, kGreen+2);

  surfCountStat.DrawContour(crit3sig, 7, kBlue);
  surfCount.DrawContour(crit3sigCoarse, kSolid, kBlue);

  TH1F* hDummy1 = new TH1F("","",1,0,1);
  TH1F* hDummy2 = new TH1F("","",1,0,1);
  TH1F* hDummy3 = new TH1F("","",1,0,1);

  hDummy1->SetLineColor(kRed);
  hDummy2->SetLineColor(kBlue);
  hDummy3->SetLineColor(kGreen+2);

  TLegend* leg = new TLegend(0.16, 0.16, 0.6/*0.4*/, .4);
  leg->SetFillStyle(0);
  leg->SetHeader("3#sigma C.L.","C");
  leg->AddEntry(hDummy2, "SBND counting only", "l");
  leg->AddEntry(hDummy3, "Icarus/SBND ratio only", "l");
  leg->AddEntry(hDummy1, "Combined fit", "l");
  //  leg->SetLineColor(10);
  //  leg->SetFillColor(10);
  leg->Draw();

  gPad->Print("nus_extrap.pdf");
}
