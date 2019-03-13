#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Loaders.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Prediction/PredictionInterp.h"
#include "CAFAna/Analysis/Calcs.h"
#include "CAFAna/Analysis/Plots.h"
#include "CAFAna/Analysis/Style.h"
#include "CAFAna/Analysis/Surface.h"
#include "CAFAna/Vars/FitVars.h"
#include "CAFAna/Experiment/SingleSampleExperiment.h"
#include "CAFAna/Experiment/MultiExperiment.h"
#include "CAFAna/Experiment/SolarConstraints.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Systs/Systs.h"
#include "CAFAna/Analysis/Fit.h"
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

void kirk_dune()
{
  //  this binning is a primary determinant of the time this macro takes to run

  //  very coarse binning, just to test macro, will worsen results: ~ 45 min
  //  int binsnumuE = 4;
  //  int binsnumuPID = 4;
  //  int binsnueE = 4;
  //  int binsnuePID = 4;

  // full binning: 2-3 hours
  int binsnumuE = 15;
  int binsnumuPID = 15;
  int binsnueE = 15;
  int binsnuePID = 20;

  int binsnumu2d = binsnumuE*binsnumuPID;
  int binsnue2d = binsnumuE*binsnumuPID;

  rootlogon(); // style

  // POT/yr * 3.5yrs * mass correction
  const double pot = 3.5 * 1.47e21 * 40/1.13;

  SpectrumLoader loaderNumu("/pnfs/dune/persistent/TaskForce_AnaTree/far/train/v2.2/numutest.root");
  SpectrumLoader loaderNue("/pnfs/dune/persistent/TaskForce_AnaTree/far/train/v2.2/nuetest.root");

  SpectrumLoader loaderNumuRHC("/pnfs/dune/persistent/TaskForce_AnaTree/far/train/v2.2/anumutest.root");
  SpectrumLoader loaderNueRHC("/pnfs/dune/persistent/TaskForce_AnaTree/far/train/v2.2/anuetest.root");

  osc::IOscCalculatorAdjustable* calc = DefaultOscCalc();
  calc->SetL(1300);
  calc->SetdCP(TMath::Pi()*1.5);

  // Standard DUNE numbers from Matt Bass
  calc->SetTh12(0.5857);
  calc->SetTh13(0.148);
  calc->SetTh23(0.7854);
  calc->SetDmsq21(0.000075);
  calc->SetDmsq32(0.002524-0.000075); // quoted value is 31

  osc::IOscCalculatorAdjustable* calci = DefaultOscCalc();
  calci->SetL(1300);
  calci->SetdCP(TMath::Pi()*1.5);

  // Standard DUNE numbers from Matt Bass
  calci->SetTh12(0.5857);
  calci->SetTh13(0.148);
  calci->SetTh23(0.7854);
  calci->SetDmsq21(0.000075);
  calci->SetDmsq32(-(0.002524+0.000075)); // quoted value is 31

  // One sigma errors
  // (t12,t13,t23,dm21,dm32)=(0.023,0.018,0.058,0.0,0.024,0.016)

  auto* loaderNumuBeam  = loaderNumu.LoaderForRunPOT(20000001);
  auto* loaderNumuNue   = loaderNumu.LoaderForRunPOT(20000002);
  auto* loaderNumuNuTau = loaderNumu.LoaderForRunPOT(20000003);
  auto* loaderNumuNC    = loaderNumu.LoaderForRunPOT(0);

  auto* loaderNueBeam  = loaderNue.LoaderForRunPOT(20000001);
  auto* loaderNueNue   = loaderNue.LoaderForRunPOT(20000002);
  auto* loaderNueNuTau = loaderNue.LoaderForRunPOT(20000003);
  auto* loaderNueNC    = loaderNue.LoaderForRunPOT(0);

  auto* loaderNumuBeamRHC  = loaderNumuRHC.LoaderForRunPOT(20000004);
  auto* loaderNumuNueRHC   = loaderNumuRHC.LoaderForRunPOT(20000005);
  auto* loaderNumuNuTauRHC = loaderNumuRHC.LoaderForRunPOT(20000006);
  auto* loaderNumuNCRHC    = loaderNumuRHC.LoaderForRunPOT(0);

  auto* loaderNueBeamRHC  = loaderNueRHC.LoaderForRunPOT(20000004);
  auto* loaderNueNueRHC   = loaderNueRHC.LoaderForRunPOT(20000005);
  auto* loaderNueNuTauRHC = loaderNueRHC.LoaderForRunPOT(20000006);
  auto* loaderNueNCRHC    = loaderNueRHC.LoaderForRunPOT(0);

  Loaders loadersdunenue;
  loadersdunenue.AddLoader(loaderNueBeam,  caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNonSwap);
  loadersdunenue.AddLoader(loaderNueNue,   caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNue);
  loadersdunenue.AddLoader(loaderNueNuTau, caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kTau);
  loadersdunenue.AddLoader(loaderNueNC,   caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNC);

  Loaders loadersdunenuerhc;
  loadersdunenuerhc.AddLoader(loaderNueBeamRHC,  caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNonSwap);
  loadersdunenuerhc.AddLoader(loaderNueNueRHC,   caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNue);
  loadersdunenuerhc.AddLoader(loaderNueNuTauRHC, caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kTau);
  loadersdunenuerhc.AddLoader(loaderNueNCRHC,   caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNC);

  Loaders loadersdunenumu;
  loadersdunenumu.AddLoader(loaderNumuBeam,  caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNonSwap);
  loadersdunenumu.AddLoader(loaderNumuNue,   caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNue);
  loadersdunenumu.AddLoader(loaderNumuNuTau, caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kTau);
  loadersdunenumu.AddLoader(loaderNumuNC,   caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNC);

  Loaders loadersdunenumurhc;
  loadersdunenumurhc.AddLoader(loaderNumuBeamRHC,  caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNonSwap);
  loadersdunenumurhc.AddLoader(loaderNumuNueRHC,   caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNue);
  loadersdunenumurhc.AddLoader(loaderNumuNuTauRHC, caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kTau);
  loadersdunenumurhc.AddLoader(loaderNumuNCRHC,   caf::kFARDET, Loaders::kMC, ana::kBeam, Loaders::kNC);

  const Var Enu_reco = SIMPLEVAR(dune.Ev_reco);
  const Var pid_reco = SIMPLEVAR(dune.mvaresult);

  float kNumuMVACutFHC = 0.25;
  float kNumuMVACutRHC = 0.5;
  float kNueMVACut = 0.8;

  const Cut kSelNumu = SIMPLEVAR(dune.mvaresult) > kNumuMVACutFHC;
  const Cut kSelNumuRHC = SIMPLEVAR(dune.mvaresult) > kNumuMVACutRHC;
  const Cut kSelNue = SIMPLEVAR(dune.mvaresult) > kNueMVACut;
  const Cut kSelNueRHC = SIMPLEVAR(dune.mvaresult) > kNueMVACut;


  std::vector<const ISyst*> systsE = {&kEnergyScaleSyst, &kEnergyResSyst};
  std::vector<const ISyst*> systsnorm = {&kNCSyst, &kNutauSyst,  &kNueBeamSyst};
  std::vector<const ISyst*> systsall = {&kEnergyScaleSyst, &kEnergyResSyst, &kNCSyst, &kNutauSyst, &kNueBeamSyst};
  std::vector<const ISyst*> systsall2 = {&kEnergyScaleSyst, &kEnergyResSyst, &kNCSyst2, &kNutauSyst, &kNueBeamSyst};


  const Var kEnuPidNumu = Var2D(Enu_reco, Binning::Simple(binsnumuE, 0, 10), pid_reco, Binning::Simple(binsnumuPID, -1, +1));
  const Var kEnuPidNue = Var2D(Enu_reco, Binning::Simple(binsnueE, 0, 6), pid_reco, Binning::Simple(binsnuePID, -1, +1));

  // 2D w/ Loaders
  NoExtrapGenerator gendunenumu(HistAxis("PID/E2d", Binning::Simple(binsnumu2d,0,binsnumu2d), kEnuPidNumu), kNoCut);
  NoExtrapGenerator gendunenue (HistAxis("PID/E2d", Binning::Simple(binsnue2d,0,binsnue2d), kEnuPidNue),  kNoCut);
  PredictionInterp preddunenumu(systsall, calc, gendunenumu, loadersdunenumu);
  PredictionInterp preddunenue (systsall, calc, gendunenue,  loadersdunenue);

  NoExtrapGenerator gendunenumurhc(HistAxis("PID/E2d", Binning::Simple(binsnumu2d,0,binsnumu2d), kEnuPidNumu), kNoCut);
  NoExtrapGenerator gendunenuerhc (HistAxis("PID/E2d", Binning::Simple(binsnue2d,0,binsnue2d), kEnuPidNue),  kNoCut);
  PredictionInterp preddunenumurhc(systsall, calc, gendunenumurhc, loadersdunenumurhc);
  PredictionInterp preddunenuerhc (systsall, calc, gendunenuerhc,  loadersdunenuerhc);


  // 1D w/ Loaders
  NoExtrapGenerator gendunenumu1d(HistAxis("Reconstructed Energy (GeV)", Binning::Simple(80,0,10), Enu_reco), kSelNumu);
  NoExtrapGenerator gendunenue1d (HistAxis("Reconstructed Energy (GeV)", Binning::Simple(80,0,10), Enu_reco), kSelNue);
  PredictionInterp preddunenumu1d(systsall, calc, gendunenumu1d, loadersdunenumu);
  PredictionInterp preddunenue1d (systsall, calc, gendunenue1d,  loadersdunenue);

  NoExtrapGenerator gendunenumu1drhc(HistAxis("Reconstructed Energy (GeV)", Binning::Simple(80,0,10), Enu_reco), kSelNumuRHC);
  NoExtrapGenerator gendunenue1drhc (HistAxis("Reconstructed Energy (GeV)", Binning::Simple(80,0,10), Enu_reco), kSelNueRHC);
  PredictionInterp preddunenumu1drhc(systsall, calc, gendunenumu1drhc, loadersdunenumurhc);
  PredictionInterp preddunenue1drhc (systsall, calc, gendunenue1drhc,  loadersdunenuerhc);


  // SpectrumLoader instead of Loaders
  PredictionNoExtrap predNumuPID(*loaderNumuBeam, *loaderNumuNue, *loaderNumuNuTau, *loaderNumuNC, "PID", Binning::Simple(100, -1, +1), SIMPLEVAR(dune.mvaresult), kNoCut);

  PredictionNoExtrap predNumuPIDRHC(*loaderNumuBeamRHC, *loaderNumuNueRHC, *loaderNumuNuTauRHC, *loaderNumuNCRHC, "PID", Binning::Simple(100, -1, +1), SIMPLEVAR(dune.mvaresult), kNoCut);

  PredictionNoExtrap predNuePID(*loaderNueBeam, *loaderNueNue, *loaderNueNuTau, *loaderNueNC, "PID", Binning::Simple(100, -1, +1), SIMPLEVAR(dune.mvaresult), kNoCut);

  PredictionNoExtrap predNuePIDRHC(*loaderNueBeamRHC, *loaderNueNueRHC, *loaderNueNuTauRHC, *loaderNueNCRHC, "PID", Binning::Simple(100, -1, +1), SIMPLEVAR(dune.mvaresult), kNoCut);

  PredictionNoExtrap pred(*loaderNumuBeam, *loaderNumuNue, *loaderNumuNuTau, *loaderNumuNC, "Reconstructed E (GeV)", Binning::Simple(80, 0, 10), Enu_reco, kSelNumu);

  PredictionNoExtrap pred2d(*loaderNumuBeam, *loaderNumuNue, *loaderNumuNuTau, *loaderNumuNC, "Reconstructed E (GeV)", Binning::Simple(binsnumu2d, 0, binsnumu2d), kEnuPidNumu, kNoCut);

  PredictionNoExtrap predRHC(*loaderNumuBeamRHC, *loaderNumuNueRHC, *loaderNumuNuTauRHC, *loaderNumuNCRHC, "Reconstructed E (GeV)", Binning::Simple(80, 0, 10), Enu_reco, kSelNumuRHC);

  PredictionNoExtrap predNue(*loaderNueBeam, *loaderNueNue, *loaderNueNuTau, *loaderNueNC, "Reconstructed E (GeV)", Binning::Simple(24, 0, 6), Enu_reco, kSelNue);

  PredictionNoExtrap predNue2d(*loaderNueBeam, *loaderNueNue, *loaderNueNuTau, *loaderNueNC, "Reconstructed E (GeV)", Binning::Simple(binsnue2d, 0, binsnue2d), kEnuPidNue, kNoCut);

  PredictionNoExtrap predNueRHC(*loaderNueBeamRHC, *loaderNueNueRHC, *loaderNueNuTauRHC, *loaderNueNCRHC, "Reconstructed E (GeV)", Binning::Simple(24, 0, 6), Enu_reco, kSelNueRHC);


  // test systematics are really shifting
  SystShifts scaleshift, resshift;
  scaleshift.SetShift(&kEnergyScaleSyst, +3);
  resshift.SetShift(&kEnergyResSyst, +3);

  Spectrum nom2(*loaderNumuBeam, HistAxis("blah", Binning::Simple(binsnumuE,0,binsnumuE), Enu_reco), kNoCut);
  Spectrum scale2(*loaderNumuBeam, HistAxis("blah", Binning::Simple(binsnumuE,0,binsnumuE), Enu_reco), kNoCut, scaleshift);
  Spectrum res2(*loaderNumuBeam, HistAxis("blah", Binning::Simple(binsnumuE,0,binsnumuE), Enu_reco), kNoCut, resshift);


  loaderNumu.Go();
  loaderNue.Go();
  loaderNumuRHC.Go();
  loaderNueRHC.Go();

  // have to have this for each prediction from a Loaders or it doesn't work; TODO fix this in source code
  preddunenue.LoadedCallback();
  preddunenumu.LoadedCallback();
  preddunenuerhc.LoadedCallback();
  preddunenumurhc.LoadedCallback();

  preddunenue1d.LoadedCallback();
  preddunenumu1d.LoadedCallback();
  preddunenue1drhc.LoadedCallback();
  preddunenumu1drhc.LoadedCallback();


  SaveToFile(preddunenue, "pred_nue_fhc.root", "pred");
  SaveToFile(preddunenuerhc, "pred_nue_rhc.root", "pred");
  SaveToFile(preddunenumu, "pred_numu_fhc.root", "pred");
  SaveToFile(preddunenumurhc, "pred_numu_rhc.root", "pred");


  Spectrum mock = pred.Predict(calc).FakeData(pot);
  SingleSampleExperiment expt(&pred, mock);

  Spectrum mock2d = pred2d.Predict(calc).FakeData(pot);
  SingleSampleExperiment expt2d(&pred2d, mock2d);

  Spectrum mock2dsysts = preddunenumu.Predict(calc).FakeData(pot);
  SingleSampleExperiment expt2dsysts(&preddunenumu, mock2dsysts);

  Spectrum mock2dsystsrhc = preddunenumurhc.Predict(calc).FakeData(pot);
  SingleSampleExperiment expt2dsystsrhc(&preddunenumurhc, mock2dsystsrhc);

  Spectrum mock2di = pred2d.Predict(calci).FakeData(pot);
  SingleSampleExperiment expt2di(&pred2d, mock2di);

  Spectrum mock2disysts = preddunenumu.Predict(calci).FakeData(pot);
  SingleSampleExperiment expt2disysts(&preddunenumu, mock2disysts);

  Spectrum mock2disystsrhc = preddunenumurhc.Predict(calci).FakeData(pot);
  SingleSampleExperiment expt2disystsrhc(&preddunenumurhc, mock2disystsrhc);

  Spectrum mockRHC = predRHC.Predict(calc).FakeData(pot);
  SingleSampleExperiment exptRHC(&predRHC, mockRHC);

  Spectrum mockNue = predNue.Predict(calc).FakeData(pot);
  SingleSampleExperiment exptNue(&predNue, mockNue);

  Spectrum mockNue2d = predNue2d.Predict(calc).FakeData(pot);
  SingleSampleExperiment exptNue2d(&predNue2d, mockNue2d);

  Spectrum mockNue2dsysts = preddunenue.Predict(calc).FakeData(pot);
  SingleSampleExperiment exptNue2dsysts(&preddunenue, mockNue2dsysts);

  Spectrum mockNue2dsystsrhc = preddunenuerhc.Predict(calc).FakeData(pot);
  SingleSampleExperiment exptNue2dsystsrhc(&preddunenuerhc, mockNue2dsystsrhc);

  Spectrum mockNue2di = predNue2d.Predict(calci).FakeData(pot);
  SingleSampleExperiment exptNue2di(&predNue2d, mockNue2di);

  Spectrum mockNue2disysts = preddunenue.Predict(calci).FakeData(pot);
  SingleSampleExperiment exptNue2disysts(&preddunenue, mockNue2disysts);

  Spectrum mockNue2disystsrhc = preddunenuerhc.Predict(calci).FakeData(pot);
  SingleSampleExperiment exptNue2disystsrhc(&preddunenuerhc, mockNue2disystsrhc);

  Spectrum mockNueRHC = predNueRHC.Predict(calc).FakeData(pot);
  SingleSampleExperiment exptNueRHC(&predNueRHC, mockNueRHC);

  Spectrum mockNuePID = predNuePID.Predict(calc).FakeData(pot);
  Spectrum mockNuePIDRHC = predNuePIDRHC.Predict(calc).FakeData(pot);
  Spectrum mockNumuPID = predNumuPID.Predict(calc).FakeData(pot);
  Spectrum mockNumuPIDRHC = predNumuPIDRHC.Predict(calc).FakeData(pot);


  Surface surf(&expt, calc, &kFitSinSqTheta23, 20, .4, .6, &kFitDmSq32Scaled, 20, 2.35, 2.55);
  Surface surf2d(&expt2d, calc, &kFitSinSqTheta23, 20, .4, .6, &kFitDmSq32Scaled, 20, 2.35, 2.55);
  Surface surf2dcheck(&expt2dsysts, calc, &kFitSinSqTheta23, 20, .4, .6, &kFitDmSq32Scaled, 20, 2.35, 2.55); //should be same as previous line, let's check
  Surface surf2dsysts(&expt2dsysts, calc, &kFitSinSqTheta23, 20, .4, .6, &kFitDmSq32Scaled, 20, 2.35, 2.55, {}, systsnorm);
  Surface surf2dsystsall(&expt2dsysts, calc, &kFitSinSqTheta23, 20, .4, .6, &kFitDmSq32Scaled, 20, 2.35, 2.55, {}, systsall);
  Surface surfRHC(&exptRHC, calc, &kFitSinSqTheta23, 20, .4, .6, &kFitDmSq32Scaled, 20, 2.35, 2.55);

  MultiExperiment numuall({&expt2dsysts, &expt2dsystsrhc});

  Surface surfAll(&numuall, calc, &kFitSinSqTheta23, 20, .4, .6, &kFitDmSq32Scaled, 20, 2.35, 2.55, {}, systsall);

  MultiExperiment me({&expt, &exptNue, new SolarConstraints()});
  MultiExperiment me2d({&expt2d, &exptNue2d, new SolarConstraints()});
  MultiExperiment me2dsysts({&expt2dsysts, &exptNue2dsysts, new SolarConstraints()});
  MultiExperiment me2di({&expt2di, &exptNue2di, new SolarConstraints()});
  MultiExperiment me2disysts({&expt2disysts, &exptNue2disysts, new SolarConstraints()});
  MultiExperiment meRHC({&exptRHC, &exptNueRHC, new SolarConstraints()});
  MultiExperiment meAll({&exptNue2dsysts, &exptNue2dsystsrhc, &expt2dsysts, &expt2dsystsrhc, new SolarConstraints()});

  Surface surfNue(&me, calc, &kFitDeltaInPiUnits, 20, 0, 2, &kFitSinSqTheta23, 20, .4, .6);
  Surface surfNue2d(&me2d, calc, &kFitDeltaInPiUnits, 20, 0, 2, &kFitSinSqTheta23, 20, .4, .6);
  Surface surfNue2dsysts(&me2dsysts, calc, &kFitDeltaInPiUnits, 20, 0, 2, &kFitSinSqTheta23, 20, .4, .6, {}, systsnorm);
  Surface surfNue2dsystsall(&me2dsysts, calc, &kFitDeltaInPiUnits, 20, 0, 2, &kFitSinSqTheta23, 20, .4, .6, {}, systsall);
  Surface surfNue2dsystsall2(&me2dsysts, calc, &kFitDeltaInPiUnits, 20, 0, 2, &kFitSinSqTheta23, 20, .4, .6, {}, systsall2);
  Surface surfNueRHC(&meRHC, calc, &kFitDeltaInPiUnits, 20, 0, 2, &kFitSinSqTheta23, 20, .4, .6);
  Surface surfNueAll(&meAll, calc, &kFitDeltaInPiUnits, 20, 0, 2, &kFitSinSqTheta23, 20, .4, .6, {}, systsall);

  new TCanvas;
  TH1* h3 = DataMCComparisonComponents(mock, &pred, calc);
  h3->SetTitle("#nu_{#mu} FHC selection (MVA>0.25) 3.5yrs #times 40kt");
  CenterTitles(h3);
  Legend();
  gPad->Print("components.pdf");
  gPad->Print("components.C");

  new TCanvas;
  TH1* h3r = DataMCComparisonComponents(mockRHC, &predRHC, calc);
  h3r->SetTitle("#nu_{#mu} RHC selection (MVA>0.5) 3.5yrs #times 40kt");
  CenterTitles(h3);
  Legend();
  gPad->Print("components_rhc.pdf");
  gPad->Print("components_rhc.C");

  new TCanvas;
  TH1* h4 = DataMCComparisonComponents(mockNue, &predNue, calc);
  h4->SetTitle("#nu_{e} FHC selection (MVA>0.8) 3.5yrs #times 40kt");
  CenterTitles(h4);
  Legend();
  gPad->Print("components_nue.pdf");
  gPad->Print("components_nue.C");

  new TCanvas;
  TH1* h4r = DataMCComparisonComponents(mockNueRHC, &predNueRHC, calc);
  h4r->SetTitle("#nu_{e} RHC selection (MVA>0.8) 3.5yrs #times 40kt");
  CenterTitles(h4);
  Legend();
  gPad->Print("components_nue_rhc.pdf");
  gPad->Print("components_nue_rhc.C");

  new TCanvas;
  TH1* h2 = DataMCComparisonComponents(mockNumuPID, &predNumuPID, calc);
  h2->SetTitle("#nu_{#mu} FHC 3.5yrs #times 40kt");
  CenterTitles(h2);
  Legend();
  h2->GetYaxis()->SetRangeUser(1, 1e4);
  gPad->SetLogy();
  gPad->Print("components_pid.pdf");
  gPad->Print("components_pid.C");

  new TCanvas;
  TH1* h2d = DataMCComparisonComponents(mock2d, &pred2d, calc);
  h2d->SetTitle("#nu_{#mu} FHC 3.5yrs #times 40kt");
  CenterTitles(h2);
  Legend();
  h2->GetYaxis()->SetRangeUser(1, 1e4);
  gPad->SetLogy();
  gPad->Print("components_pid_2d.pdf");
  gPad->Print("components_pid_2d.C");

  new TCanvas;
  TH1* h2r = DataMCComparisonComponents(mockNumuPIDRHC, &predNumuPIDRHC, calc);
  h2r->SetTitle("#nu_{#mu} RHC 3.5yrs #times 40kt");
  CenterTitles(h2);
  Legend();
  h2->GetYaxis()->SetRangeUser(1, 1e4);
  gPad->SetLogy();
  gPad->Print("components_pid_rhc.pdf");
  gPad->Print("components_pid_rhc.C");

  new TCanvas;
  TH1* h = DataMCComparisonComponents(mockNuePID, &predNuePID, calc);
  h->SetTitle("#nu_{e} FHC 3.5yrs #times 40kt");
  CenterTitles(h);
  Legend();
  h->GetYaxis()->SetRangeUser(0, 600);
  gPad->Print("components_nue_pid.pdf");
  gPad->Print("components_nue_pid.C");

  new TCanvas;
  TH1* hd = DataMCComparisonComponents(mockNue2d, &predNue2d, calc);
  hd->SetTitle("#nu_{e} FHC 3.5yrs #times 40kt");
  CenterTitles(h);
  Legend();
  h->GetYaxis()->SetRangeUser(0, 600);
  gPad->Print("components_nue_pid_2d.pdf");
  gPad->Print("components_nue_pid_2d.C");

  new TCanvas;
  TH1* hr = DataMCComparisonComponents(mockNuePIDRHC, &predNuePIDRHC, calc);
  h->SetTitle("#nu_{e} RHC 3.5yrs #times 40kt");
  CenterTitles(h);
  Legend();
  h->GetYaxis()->SetRangeUser(0, 600);
  gPad->Print("components_nue_pid_rhc.pdf");
  gPad->Print("components_nue_pid_rhc.C");

  new TCanvas;
  surf.DrawContour(Gaussian90Percent2D(surf), kSolid, 4);
  //  surfRHC.DrawContour(Gaussian90Percent2D(surfRHC), kSolid, kGreen+2);
  surf2d.DrawContour(Gaussian90Percent2D(surf2d), kSolid, 1);
  surf2dcheck.DrawContour(Gaussian90Percent2D(surf2d), kSolid, 1); // if there are 2 solid black lines there is a problem
  surf2dsysts.DrawContour(Gaussian90Percent2D(surf2dsysts), kSolid, kGreen+2);
  surf2dsystsall.DrawContour(Gaussian90Percent2D(surf2dsystsall), kSolid, 2);
  surf.DrawBestFit(kRed);
  surfAll.DrawContour(Gaussian90Percent2D(surfAll), 1, kMagenta);
  gPad->Print("cont.pdf");

  new TCanvas;
  surfNue.DrawContour(Gaussian90Percent2D(surfNue), kSolid, 4);
  surfNue2d.DrawContour(Gaussian90Percent2D(surfNue2d), kSolid, 1);
  surfNue2dsysts.DrawContour(Gaussian90Percent2D(surfNue2dsysts), kSolid, kGreen+2);
  surfNue2dsystsall.DrawContour(Gaussian90Percent2D(surfNue2dsystsall), kSolid, 2);
  //  surfNue2dsystsall2.DrawContour(Gaussian90Percent2D(surfNue2dsystsall2), kSolid, kMagenta); // just a check that a 50% NC syst is same as a 5%
  surfNueAll.DrawContour(Gaussian90Percent2D(surfNueAll), kSolid, kMagenta);
  surfNue.DrawBestFit(kRed);
  gPad->Print("cont_nue.pdf");

  Spectrum testnom = preddunenumu.Predict(calc);
  Spectrum scaleshifted = preddunenumu.PredictSyst(calc,scaleshift);
  Spectrum resshifted = preddunenumu.PredictSyst(calc,resshift);

  new TCanvas;

  TH1* htestnom = testnom.ToTH1(pot);
  TH1* hscaleshifted = scaleshifted.ToTH1(pot);
  TH1* hresshifted = resshifted.ToTH1(pot);

  htestnom->SetLineWidth(2);
  hscaleshifted->SetLineWidth(2);
  hresshifted->SetLineWidth(2);

  hscaleshifted->SetLineColor(2);
  hresshifted->SetLineColor(4);

  htestnom->Draw();
  hscaleshifted->Draw("same");
  hresshifted->Draw("same");

  gPad->Print("testsysts.pdf"); // yes things are really shifting, here we can see it happen

  new TCanvas;

  TH1* hnom2 = nom2.ToTH1(pot);
  TH1* hscale2 = scale2.ToTH1(pot);
  TH1* hres2 = res2.ToTH1(pot);

  hnom2->SetLineWidth(2);
  hscale2->SetLineWidth(2);
  hres2->SetLineWidth(2);

  hscale2->SetLineColor(2);
  hres2->SetLineColor(4);

  hnom2->Draw();
  hscale2->Draw("same");
  hres2->Draw("same");

  gPad->Print("testsysts2.pdf"); // check shifts are happening in 2D variable also

  new TCanvas;

  // This is a very cheesy way to make the McD plot - would have to be very
  // different if we were varying any other parameters; leave in for now but don't trust this
  calc->SetdCP(0);
  Spectrum zeroNumu = pred.Predict(calc).FakeData(pot);
  Spectrum zeroNue = predNue.Predict(calc).FakeData(pot);
  Spectrum zeroNumu2d = pred2d.Predict(calc).FakeData(pot);
  Spectrum zeroNue2d = predNue2d.Predict(calc).FakeData(pot);
  Spectrum zeroNumuRHC = predRHC.Predict(calc).FakeData(pot);
  Spectrum zeroNueRHC = predNueRHC.Predict(calc).FakeData(pot);
  calc->SetdCP(TMath::Pi());
  Spectrum oneNumu = pred.Predict(calc).FakeData(pot);
  Spectrum oneNue = predNue.Predict(calc).FakeData(pot);
  Spectrum oneNumu2d = pred2d.Predict(calc).FakeData(pot);
  Spectrum oneNue2d = predNue2d.Predict(calc).FakeData(pot);
  Spectrum oneNumuRHC = predRHC.Predict(calc).FakeData(pot);
  Spectrum oneNueRHC = predNueRHC.Predict(calc).FakeData(pot);

  TGraph* g = new TGraph;
  TGraph* g2d = new TGraph;
  TGraph* grhc = new TGraph;
  TGraph* gold = new TGraph;
  TGraph* gold2 = new TGraph;

  for(int i = -100; i <= 100; ++i){
    calc->SetdCP(i/100.*TMath::Pi());

    Spectrum mockNumu = pred.Predict(calc).FakeData(pot);
    Spectrum mockNue = predNue.Predict(calc).FakeData(pot);
    Spectrum mockNumu2d = pred2d.Predict(calc).FakeData(pot);
    Spectrum mockNue2d = predNue2d.Predict(calc).FakeData(pot);
    Spectrum mockNumuRHC = predRHC.Predict(calc).FakeData(pot);
    Spectrum mockNueRHC = predNueRHC.Predict(calc).FakeData(pot);

    const double llZero = LogLikelihood(zeroNumu.ToTH1(pot), mockNumu.ToTH1(pot))+
      LogLikelihood(zeroNue.ToTH1(pot), mockNue.ToTH1(pot));

    const double llOne = LogLikelihood(oneNumu.ToTH1(pot), mockNumu.ToTH1(pot))+
      LogLikelihood(oneNue.ToTH1(pot), mockNue.ToTH1(pot));

    const double llZero2d = LogLikelihood(zeroNumu2d.ToTH1(pot), mockNumu2d.ToTH1(pot))+
      LogLikelihood(zeroNue2d.ToTH1(pot), mockNue2d.ToTH1(pot));

    const double llOne2d = LogLikelihood(oneNumu2d.ToTH1(pot), mockNumu2d.ToTH1(pot))+
      LogLikelihood(oneNue2d.ToTH1(pot), mockNue2d.ToTH1(pot));

    const double llZeroRHC = LogLikelihood(zeroNumuRHC.ToTH1(pot), mockNumuRHC.ToTH1(pot))+LogLikelihood(zeroNumu.ToTH1(pot), mockNumu.ToTH1(pot))+
      LogLikelihood(zeroNueRHC.ToTH1(pot), mockNueRHC.ToTH1(pot))+LogLikelihood(zeroNue.ToTH1(pot), mockNue.ToTH1(pot));

    const double llOneRHC = LogLikelihood(oneNumuRHC.ToTH1(pot), mockNumuRHC.ToTH1(pot))+LogLikelihood(oneNumu.ToTH1(pot), mockNumu.ToTH1(pot))+
      LogLikelihood(oneNueRHC.ToTH1(pot), mockNueRHC.ToTH1(pot))+LogLikelihood(oneNue.ToTH1(pot), mockNue.ToTH1(pot));

    const double ll = std::min(llZero, llOne);
    const double ll2d = std::min(llZero2d, llOne2d);
    const double llrhc = std::min(llZeroRHC, llOneRHC);

    if(ll>0) g->SetPoint(g->GetN(), i/100., sqrt(ll));
    if(ll2d>0) g2d->SetPoint(g2d->GetN(), i/100., sqrt(ll2d));
    if(llrhc>0) grhc->SetPoint(grhc->GetN(), i/100., sqrt(llrhc));
  }

  TH2* axes = new TH2F("", "3.5yrs #times 40kt, stats only;#delta / #pi;#sqrt{#chi^{2}}", 200, -1, +1, 50, 0, 8);

  CenterTitles(axes);
  axes->Draw();

  g2d->SetLineWidth(2);
  g2d->SetLineColor(1);
  g2d->Draw("l same");

  g->SetLineWidth(2);
  g->SetLineColor(4);
  g->Draw("l same");

  grhc->SetLineWidth(2);
  grhc->SetLineColor(kMagenta);
  grhc->Draw("l same");

  gPad->SetGridx();
  gPad->SetGridy();

  gPad->Print("mcd_hacky.pdf");
  gPad->Print("mcd_hacky.C");



  new TCanvas; // now try to make the plot right

  TH2* axes2 = new TH2F("", "3.5yrs #times 40kt, stats only;#delta / #pi;#sqrt{#chi^{2}}", 100, -1, +1, 50, 0, 8);

  TGraph* gr1d = new TGraph;    // IH+NH, all systs, 1D (E only)
  TGraph* grfull = new TGraph;  // IH+NH, all systs, 2D (ExPID)
  TGraph* gr2 = new TGraph;     // IH+NH, norm systs only, 2D
  TGraph* gr3 = new TGraph;     // FHC only, all systs, 2D, IH+NH
  TGraph* grstat = new TGraph;  // FHC only, no  systs, 2D, IH+NH
  TGraph* grnh = new TGraph;    // NH only part of previous
  TGraph* grih = new TGraph;    // IH only part


  // lots of possible combinations here: 2 hierarchies x 2 (FHC, RHC) x 2 expts (numu, nue) 2 x variables (1D, 2D) = 16, at 2 values of dCP = 32

  // do dCP = 0 case

  calc->SetdCP(0);
  calci->SetdCP(0);

  Spectrum tempnumu = preddunenumu.Predict(calc).FakeData(pot);
  SingleSampleExperiment expttempnumu(&preddunenumu, tempnumu);
  Spectrum tempnue = preddunenue.Predict(calc).FakeData(pot);
  SingleSampleExperiment expttempnue(&preddunenue, tempnue);

  Spectrum tempnumui = preddunenumu.Predict(calci).FakeData(pot);
  SingleSampleExperiment expttempnumui(&preddunenumu, tempnumu);
  Spectrum tempnuei = preddunenue.Predict(calci).FakeData(pot);
  SingleSampleExperiment expttempnuei(&preddunenue, tempnue);

  Spectrum tempnumurhc = preddunenumurhc.Predict(calc).FakeData(pot);
  SingleSampleExperiment expttempnumurhc(&preddunenumurhc, tempnumurhc);
  Spectrum tempnuerhc = preddunenuerhc.Predict(calc).FakeData(pot);
  SingleSampleExperiment expttempnuerhc(&preddunenuerhc, tempnuerhc);

  Spectrum tempnumuirhc = preddunenumurhc.Predict(calci).FakeData(pot);
  SingleSampleExperiment expttempnumuirhc(&preddunenumurhc, tempnumurhc);
  Spectrum tempnueirhc = preddunenuerhc.Predict(calci).FakeData(pot);
  SingleSampleExperiment expttempnueirhc(&preddunenuerhc, tempnuerhc);

  Spectrum tempnumu1d = preddunenumu1d.Predict(calc).FakeData(pot);
  SingleSampleExperiment expttempnumu1d(&preddunenumu1d, tempnumu1d);
  Spectrum tempnue1d = preddunenue1d.Predict(calc).FakeData(pot);
  SingleSampleExperiment expttempnue1d(&preddunenue1d, tempnue1d);

  Spectrum tempnumui1d = preddunenumu1d.Predict(calci).FakeData(pot);
  SingleSampleExperiment expttempnumui1d(&preddunenumu1d, tempnumu1d);
  Spectrum tempnuei1d = preddunenue1d.Predict(calci).FakeData(pot);
  SingleSampleExperiment expttempnuei1d(&preddunenue1d, tempnue1d);

  Spectrum tempnumurhc1d = preddunenumu1drhc.Predict(calc).FakeData(pot);
  SingleSampleExperiment expttempnumurhc1d(&preddunenumu1drhc, tempnumurhc1d);
  Spectrum tempnuerhc1d = preddunenue1drhc.Predict(calc).FakeData(pot);
  SingleSampleExperiment expttempnuerhc1d(&preddunenue1drhc, tempnuerhc1d);

  Spectrum tempnumuirhc1d = preddunenumu1drhc.Predict(calci).FakeData(pot);
  SingleSampleExperiment expttempnumuirhc1d(&preddunenumu1drhc, tempnumurhc1d);
  Spectrum tempnueirhc1d = preddunenue1drhc.Predict(calci).FakeData(pot);
  SingleSampleExperiment expttempnueirhc1d(&preddunenue1drhc, tempnuerhc1d);

  MultiExperiment tempme({&expttempnumu, &expttempnue, new SolarConstraints()}); // FHC only NH numu+nue
  MultiExperiment tempmei({&expttempnumui, &expttempnuei, new SolarConstraints()}); // FHC only IH numu+nue

  MultiExperiment tempmefull({&expttempnumu, &expttempnue, &expttempnumurhc, &expttempnuerhc, new SolarConstraints()}); // nue+numu FHC+RHC NH 2D
  MultiExperiment tempmeifull({&expttempnumui, &expttempnuei, &expttempnumuirhc, &expttempnueirhc, new SolarConstraints()}); // nue+numu FHC+RHC IH 2D

  MultiExperiment tempmefull1d({&expttempnumu1d, &expttempnue1d, &expttempnumurhc1d, &expttempnuerhc1d, new SolarConstraints()}); // nue+numu FHC+RHC NH 1D
  MultiExperiment tempmeifull1d({&expttempnumui1d, &expttempnuei1d, &expttempnumuirhc1d, &expttempnueirhc1d, new SolarConstraints()}); // nue+numu FHC+RHC IH 1D

  // now dCP = PI case

  calc->SetdCP(TMath::Pi());
  calci->SetdCP(TMath::Pi());

  Spectrum tempnumu2 = preddunenumu.Predict(calc).FakeData(pot);
  SingleSampleExperiment expttempnumu2(&preddunenumu, tempnumu2);
  Spectrum tempnue2 = preddunenue.Predict(calc).FakeData(pot);
  SingleSampleExperiment expttempnue2(&preddunenue, tempnue2);

  Spectrum tempnumui2 = preddunenumu.Predict(calci).FakeData(pot);
  SingleSampleExperiment expttempnumui2(&preddunenumu, tempnumu2);
  Spectrum tempnuei2 = preddunenue.Predict(calci).FakeData(pot);
  SingleSampleExperiment expttempnuei2(&preddunenue, tempnue2);

  Spectrum tempnumurhc2 = preddunenumurhc.Predict(calc).FakeData(pot);
  SingleSampleExperiment expttempnumurhc2(&preddunenumurhc, tempnumurhc2);
  Spectrum tempnuerhc2 = preddunenuerhc.Predict(calc).FakeData(pot);
  SingleSampleExperiment expttempnuerhc2(&preddunenuerhc, tempnuerhc2);

  Spectrum tempnumuirhc2 = preddunenumurhc.Predict(calci).FakeData(pot);
  SingleSampleExperiment expttempnumuirhc2(&preddunenumurhc, tempnumurhc2);
  Spectrum tempnueirhc2 = preddunenuerhc.Predict(calci).FakeData(pot);
  SingleSampleExperiment expttempnueirhc2(&preddunenuerhc, tempnuerhc2);

  Spectrum tempnumu1d2 = preddunenumu1d.Predict(calc).FakeData(pot);
  SingleSampleExperiment expttempnumu1d2(&preddunenumu1d, tempnumu1d2);
  Spectrum tempnue1d2 = preddunenue1d.Predict(calc).FakeData(pot);
  SingleSampleExperiment expttempnue1d2(&preddunenue1d, tempnue1d2);

  Spectrum tempnumui1d2 = preddunenumu1d.Predict(calci).FakeData(pot);
  SingleSampleExperiment expttempnumui1d2(&preddunenumu1d, tempnumu1d2);
  Spectrum tempnuei1d2 = preddunenue1d.Predict(calci).FakeData(pot);
  SingleSampleExperiment expttempnuei1d2(&preddunenue1d, tempnue1d2);

  Spectrum tempnumurhc1d2 = preddunenumu1drhc.Predict(calc).FakeData(pot);
  SingleSampleExperiment expttempnumurhc1d2(&preddunenumu1drhc, tempnumurhc1d2);
  Spectrum tempnuerhc1d2 = preddunenue1drhc.Predict(calc).FakeData(pot);
  SingleSampleExperiment expttempnuerhc1d2(&preddunenue1drhc, tempnuerhc1d2);

  Spectrum tempnumuirhc1d2 = preddunenumu1drhc.Predict(calci).FakeData(pot);
  SingleSampleExperiment expttempnumuirhc1d2(&preddunenumu1drhc, tempnumurhc1d2);
  Spectrum tempnueirhc1d2 = preddunenue1drhc.Predict(calci).FakeData(pot);
  SingleSampleExperiment expttempnueirhc1d2(&preddunenue1drhc, tempnuerhc1d2);

  MultiExperiment tempme2({&expttempnumu2, &expttempnue2, new SolarConstraints()});
  MultiExperiment tempmei2({&expttempnumui2, &expttempnuei2, new SolarConstraints()});

  MultiExperiment tempmefull2({&expttempnumu2, &expttempnue2, &expttempnumurhc2, &expttempnuerhc2, new SolarConstraints()});
  MultiExperiment tempmeifull2({&expttempnumui2, &expttempnuei2, &expttempnumuirhc2, &expttempnueirhc2, new SolarConstraints()});

  MultiExperiment tempmefull1d2({&expttempnumu1d2, &expttempnue1d2, &expttempnumurhc1d2, &expttempnuerhc1d2, new SolarConstraints()});
  MultiExperiment tempmeifull1d2({&expttempnumui1d2, &expttempnuei1d2, &expttempnumuirhc1d2, &expttempnueirhc1d2, new SolarConstraints()});

  // now make fitters that will do actual fits // TODO make variable names intelligible
  // these allow Dmsq32, sinsqth23, sinsq2th13 to float, constrained only by the data itself
  // they also allow Dmsq21 and sinsq2th12 to float, constrained be external solar data (see CAFAna/Experiments/SolarConstraints.*)

  // dCP = 0 fits
  Fitter fit(&tempme, {&kFitDmSq32Scaled, &kFitSinSq2Theta23, &kFitSinSq2Theta13, &kFitSinSq2Theta12, &kFitDmSq21}, systsnorm);            // FHC only, NH, norm systs only
  Fitter fit2(&tempme, {&kFitDmSq32Scaled, &kFitSinSq2Theta23, &kFitSinSq2Theta13, &kFitSinSq2Theta12, &kFitDmSq21});                      // FHC only, NH, no systs
  Fitter fit3(&tempme, {&kFitDmSq32Scaled, &kFitSinSq2Theta23, &kFitSinSq2Theta13, &kFitSinSq2Theta12, &kFitDmSq21}, systsall);            // FHC only, NH, all systs
  Fitter fitfull(&tempmefull, {&kFitDmSq32Scaled, &kFitSinSq2Theta23, &kFitSinSq2Theta13, &kFitSinSq2Theta12, &kFitDmSq21}, systsall);     // FHC+ RHC, NH, all systs
  Fitter fitfull1d(&tempmefull1d, {&kFitDmSq32Scaled, &kFitSinSq2Theta23, &kFitSinSq2Theta13, &kFitSinSq2Theta12, &kFitDmSq21}, systsall); // FHC+ RHC, NH, all systs, 1D
  SystShifts systsout = SystShifts::Nominal();
  Fitter fiti(&tempmei, {&kFitDmSq32Scaled, &kFitSinSq2Theta23, &kFitSinSq2Theta13, &kFitSinSq2Theta12, &kFitDmSq21}, systsnorm);          // same but IH now
  Fitter fit2i(&tempmei, {&kFitDmSq32Scaled, &kFitSinSq2Theta23, &kFitSinSq2Theta13, &kFitSinSq2Theta12, &kFitDmSq21});
  Fitter fit3i(&tempmei, {&kFitDmSq32Scaled, &kFitSinSq2Theta23, &kFitSinSq2Theta13, &kFitSinSq2Theta12, &kFitDmSq21}, systsall);
  Fitter fitfulli(&tempmeifull, {&kFitDmSq32Scaled, &kFitSinSq2Theta23, &kFitSinSq2Theta13, &kFitSinSq2Theta12, &kFitDmSq21}, systsall);
  Fitter fitfull1di(&tempmeifull1d, {&kFitDmSq32Scaled, &kFitSinSq2Theta23, &kFitSinSq2Theta13, &kFitSinSq2Theta12, &kFitDmSq21}, systsall);
  SystShifts systsouti = SystShifts::Nominal();


  // dCP = PI fits
  Fitter xfit(&tempme2, {&kFitDmSq32Scaled, &kFitSinSq2Theta23, &kFitSinSq2Theta13, &kFitSinSq2Theta12, &kFitDmSq21}, systsnorm);
  Fitter xfit2(&tempme2, {&kFitDmSq32Scaled, &kFitSinSq2Theta23, &kFitSinSq2Theta13, &kFitSinSq2Theta12, &kFitDmSq21});
  Fitter xfit3(&tempme2, {&kFitDmSq32Scaled, &kFitSinSq2Theta23, &kFitSinSq2Theta13, &kFitSinSq2Theta12, &kFitDmSq21}, systsall);
  Fitter xfitfull(&tempmefull2, {&kFitDmSq32Scaled, &kFitSinSq2Theta23, &kFitSinSq2Theta13, &kFitSinSq2Theta12, &kFitDmSq21}, systsall);
  Fitter xfitfull1d(&tempmefull1d2, {&kFitDmSq32Scaled, &kFitSinSq2Theta23, &kFitSinSq2Theta13, &kFitSinSq2Theta12, &kFitDmSq21}, systsall);
  SystShifts xsystsout = SystShifts::Nominal();
  Fitter xfiti(&tempmei2, {&kFitDmSq32Scaled, &kFitSinSq2Theta23, &kFitSinSq2Theta13, &kFitSinSq2Theta12, &kFitDmSq21}, systsnorm);
  Fitter xfit2i(&tempmei2, {&kFitDmSq32Scaled, &kFitSinSq2Theta23, &kFitSinSq2Theta13, &kFitSinSq2Theta12, &kFitDmSq21});
  Fitter xfit3i(&tempmei2, {&kFitDmSq32Scaled, &kFitSinSq2Theta23, &kFitSinSq2Theta13, &kFitSinSq2Theta12, &kFitDmSq21}, systsall);
  Fitter xfitfulli(&tempmeifull2, {&kFitDmSq32Scaled, &kFitSinSq2Theta23, &kFitSinSq2Theta13, &kFitSinSq2Theta12, &kFitDmSq21}, systsall);
  Fitter xfitfull1di(&tempmeifull1d2, {&kFitDmSq32Scaled, &kFitSinSq2Theta23, &kFitSinSq2Theta13, &kFitSinSq2Theta12, &kFitDmSq21}, systsall);
  SystShifts xsystsouti = SystShifts::Nominal();

  for(int i = -50; i <= 50; ++i){    // scan over dCP, create fake data at each dCP value, compare to dCP = 0 or PI cases, for each hierarchy
    calc->SetdCP(i/50.*TMath::Pi());
    calci->SetdCP(i/50.*TMath::Pi());

    std::cout << "Doing fits for dCP significance: " << i+50 << "%" << std::endl;

    // do the actual fits here, returns are chisq value of fits
    double throwaway = fit.Fit(calc, systsout, Fitter::kQuiet);  // debugging some odd behavior
    double chisqout = fit.Fit(calc, systsout, Fitter::kQuiet);
    double chisqouti = fiti.Fit(calci, systsouti, Fitter::kQuiet);
    double chisqout2 = fit2.Fit(calc, systsout, Fitter::kQuiet);
    double chisqout2i = fit2i.Fit(calci, systsout, Fitter::kQuiet);
    double chisqout3 = fit3.Fit(calc, systsout, Fitter::kQuiet);
    double chisqout3i = fit3i.Fit(calci, systsouti, Fitter::kQuiet);
    double chisqoutfull = fitfull.Fit(calc, systsout, Fitter::kQuiet);
    double chisqoutfulli = fitfulli.Fit(calci, systsout, Fitter::kQuiet);
    double chisqoutfull1d = fitfull1d.Fit(calc, systsout, Fitter::kQuiet);
    double chisqoutfull1di = fitfull1di.Fit(calci, systsout, Fitter::kQuiet);

    double xthrowaway = xfit.Fit(calc, systsout, Fitter::kQuiet);
    double xchisqout = xfit.Fit(calc, systsout, Fitter::kQuiet);
    double xchisqouti = xfiti.Fit(calci, systsouti, Fitter::kQuiet);
    double xchisqout2 = xfit2.Fit(calc, systsout, Fitter::kQuiet);
    double xchisqout2i = xfit2i.Fit(calci, systsout, Fitter::kQuiet);
    double xchisqout3 = xfit3.Fit(calc, systsout, Fitter::kQuiet);
    double xchisqout3i = xfit3i.Fit(calci, systsouti, Fitter::kQuiet);
    double xchisqoutfull = xfitfull.Fit(calc, systsout, Fitter::kQuiet);
    double xchisqoutfulli = xfitfulli.Fit(calci, systsout, Fitter::kQuiet);
    double xchisqoutfull1d = xfitfull1d.Fit(calc, systsout, Fitter::kQuiet);
    double xchisqoutfull1di = xfitfull1di.Fit(calci, systsout, Fitter::kQuiet);
    //    std::cout << " for dcp = " << i/50.*TMath::Pi() << " chisqs are: " << chisqout << " " << chisqouti << " " << sqrt(std::min(chisqout,chisqouti)) << std::endl; // debug

    // prevent possible errors
    if(chisqout < 0) chisqout = 0;
    if(chisqouti < 0) chisqouti = 0;
    if(chisqout2 < 0) chisqout = 0;
    if(chisqout2i < 0) chisqouti = 0;
    if(chisqout3 < 0) chisqout = 0;
    if(chisqout3i < 0) chisqouti = 0;
    if(chisqoutfull < 0) chisqoutfull = 0;
    if(chisqoutfulli < 0) chisqoutfulli = 0;
    if(chisqoutfull1d < 0) chisqoutfull1d = 0;
    if(chisqoutfull1di < 0) chisqoutfull1di = 0;

    if(xchisqout < 0) xchisqout = 0;
    if(xchisqouti < 0) xchisqouti = 0;
    if(xchisqout2 < 0) xchisqout = 0;
    if(xchisqout2i < 0) xchisqouti = 0;
    if(xchisqout3 < 0) xchisqout = 0;
    if(xchisqout3i < 0) xchisqouti = 0;
    if(xchisqoutfull < 0) xchisqoutfull = 0;
    if(xchisqoutfulli < 0) xchisqoutfulli = 0;
    if(xchisqoutfull1d < 0) xchisqoutfull1d = 0;
    if(xchisqoutfull1di < 0) xchisqoutfull1di = 0;

    //chisq for each point is minimum chisq of fit comparing to dCP = 0, PI, for either hierarchy

    gr2->SetPoint(gr2->GetN(), i/50., sqrt(std::min( std::min(chisqout,chisqouti), std::min(xchisqout,xchisqouti) )));
    gr3->SetPoint(gr3->GetN(), i/50., sqrt(std::min(std::min(chisqout3,chisqout3i),std::min(xchisqout3,xchisqout3i))));
    grnh->SetPoint(grnh->GetN(), i/50., sqrt(std::min(chisqout2,xchisqout2)));
    grih->SetPoint(grih->GetN(), i/50., sqrt(std::min(chisqout2i,xchisqout2i)));
    grstat->SetPoint(grstat->GetN(), i/50., sqrt(std::min(std::min(chisqout2i,chisqout2),std::min(xchisqout2i,xchisqout2))));
    grfull->SetPoint(grfull->GetN(), i/50., sqrt(std::min(std::min(chisqoutfulli,chisqoutfull),std::min(xchisqoutfulli,xchisqoutfull))));
    gr1d->SetPoint(gr1d->GetN(), i/50., sqrt(std::min(std::min(chisqoutfull1di,chisqoutfull1d),std::min(xchisqoutfull1di,xchisqoutfull1d))));
  }

  axes2->Draw();

  grfull->SetLineWidth(2);
  grfull->SetLineColor(2);
  grfull->Draw("l same");

  gr1d->SetLineWidth(2);
  gr1d->SetLineColor(4);
  gr1d->Draw("l same");

  gr2->SetLineWidth(2);
  gr2->SetLineColor(4);
  //  gr2->Draw("l same");

  gr3->SetLineWidth(2);
  gr3->SetLineColor(kGreen+2);
  gr3->Draw("l same");

  grnh->SetLineWidth(2);
  grnh->SetLineColor(1);
  //  grnh->Draw("l same");

  grih->SetLineWidth(2);
  grih->SetLineColor(2);
  //  grih->Draw("l same");

  grstat->SetLineWidth(2);
  grstat->SetLineColor(1);
  grstat->Draw("l same");

  gPad->SetGridx();
  gPad->SetGridy();

  gPad->Print("mcd.pdf");
  gPad->Print("mcd.C");
}
