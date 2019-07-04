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

void kirk_dune_quick()
{
  //  this binning is a primary determinant of the time this macro takes to run

  //  very coarse binning, just to test macro, will worsen results: ~ 15 min
  //  int binsnumuE = 4;
  //  int binsnumuPID = 4;
  //  int binsnueE = 4;
  //  int binsnuePID = 4;

  // full binning:  ~1 hour
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

  std::vector<const ISyst*> systs = {&kEnergyScaleSyst, &kEnergyResSyst, &kNCSyst, &kNutauSyst, &kNueBeamSyst};

  const Var kEnuPidNumu = Var2D(Enu_reco, Binning::Simple(binsnumuE, 0, 10), pid_reco, Binning::Simple(binsnumuPID, -1, +1));
  const Var kEnuPidNue = Var2D(Enu_reco, Binning::Simple(binsnueE, 0, 6), pid_reco, Binning::Simple(binsnuePID, -1, +1));

  // 2D w/ Loaders
  NoExtrapGenerator gendunenumu(HistAxis("PID/E2d", Binning::Simple(binsnumu2d,0,binsnumu2d), kEnuPidNumu), kNoCut);
  NoExtrapGenerator gendunenue (HistAxis("PID/E2d", Binning::Simple(binsnue2d,0,binsnue2d), kEnuPidNue),  kNoCut);
  PredictionInterp preddunenumu(systs, calc, gendunenumu, loadersdunenumu);
  PredictionInterp preddunenue (systs, calc, gendunenue,  loadersdunenue);

  NoExtrapGenerator gendunenumurhc(HistAxis("PID/E2d", Binning::Simple(binsnumu2d,0,binsnumu2d), kEnuPidNumu), kNoCut);
  NoExtrapGenerator gendunenuerhc (HistAxis("PID/E2d", Binning::Simple(binsnue2d,0,binsnue2d), kEnuPidNue),  kNoCut);
  PredictionInterp preddunenumurhc(systs, calc, gendunenumurhc, loadersdunenumurhc);
  PredictionInterp preddunenuerhc (systs, calc, gendunenuerhc,  loadersdunenuerhc);

  loaderNumu.Go();
  loaderNue.Go();
  loaderNumuRHC.Go();
  loaderNueRHC.Go();

  // have to have this for each prediction from a Loaders or it doesn't work; TODO fix this in source code
  preddunenue.LoadedCallback();
  preddunenumu.LoadedCallback();
  preddunenuerhc.LoadedCallback();
  preddunenumurhc.LoadedCallback();

  new TCanvas;

  TH2* axes2 = new TH2F("", "3.5yrs #times 40kt, stats only;#delta / #pi;#sqrt{#chi^{2}}", 100, -1, +1, 50, 0, 8);

  TGraph* grfull = new TGraph;  // IH+NH, all systs, 2D (ExPID)

  // lots of possible combinations here: 2 hierarchies x 2 (FHC, RHC) x 2 expts (numu, nue) 2 = 8, at 2 values of dCP = 16

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

  MultiExperiment tempme({&expttempnumu, &expttempnue, &expttempnumurhc, &expttempnuerhc, new SolarConstraints()}); // nue+numu FHC+RHC NH 2D
  MultiExperiment tempmei({&expttempnumui, &expttempnuei, &expttempnumuirhc, &expttempnueirhc, new SolarConstraints()}); // nue+numu FHC+RHC IH 2D

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

  MultiExperiment tempme2({&expttempnumu2, &expttempnue2, &expttempnumurhc2, &expttempnuerhc2, new SolarConstraints()});
  MultiExperiment tempmei2({&expttempnumui2, &expttempnuei2, &expttempnumuirhc2, &expttempnueirhc2, new SolarConstraints()});

  // now make fitters that will do actual fits // TODO make variable names intelligible
  // these allow Dmsq32, sinsqth23, sinsq2th13 to float, constrained only by the data itself
  // they also allow Dmsq21 and sinsq2th12 to float, constrained be external solar data (see CAFAna/Experiments/SolarConstraints.*)

  // dCP = 0 fits
  Fitter fit(&tempme, {&kFitDmSq32Scaled, &kFitSinSq2Theta23, &kFitSinSq2Theta13, &kFitSinSq2Theta12, &kFitDmSq21}, systs);     // FHC+ RHC, NH, all systs
  SystShifts systsout = SystShifts::Nominal();
  Fitter fiti(&tempmei, {&kFitDmSq32Scaled, &kFitSinSq2Theta23, &kFitSinSq2Theta13, &kFitSinSq2Theta12, &kFitDmSq21}, systs);
  SystShifts systsouti = SystShifts::Nominal();

  // dCP = PI fits
  Fitter xfit(&tempme2, {&kFitDmSq32Scaled, &kFitSinSq2Theta23, &kFitSinSq2Theta13, &kFitSinSq2Theta12, &kFitDmSq21}, systs);
  SystShifts xsystsout = SystShifts::Nominal();
  Fitter xfiti(&tempmei2, {&kFitDmSq32Scaled, &kFitSinSq2Theta23, &kFitSinSq2Theta13, &kFitSinSq2Theta12, &kFitDmSq21}, systs);
  SystShifts xsystsouti = SystShifts::Nominal();

  for(int i = -50; i <= 50; ++i){    // scan over dCP, create fake data at each dCP value, compare to dCP = 0 or PI cases, for each hierarchy
    calc->SetdCP(i/50.*TMath::Pi());
    calci->SetdCP(i/50.*TMath::Pi());

    std::cout << "Doing fits for dCP significance: " << i+50 << "%" << std::endl;

    // do the actual fits here, returns are chisq value of fits
    //    double throwaway = fit.Fit(calc, systsout, Fitter::kQuiet);  // debugging some odd behavior
    double chisqout = fit.Fit(calc, systsout, Fitter::kQuiet);
    double chisqouti = fiti.Fit(calci, systsout, Fitter::kQuiet);

    //    double xthrowaway = xfit.Fit(calc, systsout, Fitter::kQuiet);
    double xchisqout = xfit.Fit(calc, systsout, Fitter::kQuiet);
    double xchisqouti = xfiti.Fit(calci, systsout, Fitter::kQuiet);
    //    std::cout << " for dcp = " << i/50.*TMath::Pi() << " chisqs are: " << chisqout << " " << chisqouti << " " << sqrt(std::min(chisqout,chisqouti)) << std::endl; // debug

    // prevent possible errors
    if(chisqout < 0) chisqout = 0;
    if(chisqouti < 0) chisqouti = 0;
    if(xchisqout < 0) xchisqout = 0;
    if(xchisqouti < 0) xchisqouti = 0;

    //chisq for each point is minimum chisq of fit comparing to dCP = 0, PI, for either hierarchy

    grfull->SetPoint(grfull->GetN(), i/50., sqrt(std::min(std::min(chisqouti,chisqout),std::min(xchisqouti,xchisqout))));
  }

  axes2->Draw();

  grfull->SetLineWidth(2);
  grfull->SetLineColor(1);
  grfull->Draw("l same");

  gPad->SetGridx();
  gPad->SetGridy();

  gPad->Print("mcd_quick.pdf");
  gPad->Print("mcd_quick.C");
}
