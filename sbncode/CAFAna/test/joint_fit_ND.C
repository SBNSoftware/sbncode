#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Cuts/TruthCuts.h"

#include "CAFAna/Systs/DUNEXSecSysts.h"
#include "CAFAna/Systs/DUNEFluxSysts.h"
#include "CAFAna/Systs/DUNENDSysts.h"
#include "CAFAna/Core/Loaders.h"
#include "CAFAna/Core/Ratio.h"
#include "CAFAna/Core/SystShifts.h"
#include "CAFAna/Experiment/SingleSampleExperiment.h"
#include "CAFAna/Experiment/MultiExperiment.h"
#include "CAFAna/Analysis/Calcs.h"
#include "CAFAna/Analysis/Fit.h"
#include "CAFAna/Analysis/Plots.h"
#include "CAFAna/Vars/FitVars.h"
#include "CAFAna/Prediction/PredictionInterp.h"
#include "CAFAna/Prediction/PredictionNoOsc.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"


using namespace ana;

#include "Utilities/rootlogon.C"

#include "StandardRecord/StandardRecord.h"

#include "TCanvas.h"
#include "TH1.h"
#include "TGraph.h"
#include "TLegend.h"
#include "THStack.h"

const Var kRecoE = SIMPLEVAR(dune.Ev_reco);
const Var kPIDmu = SIMPLEVAR(dune.numu_pid);
const Var kPIDe = SIMPLEVAR(dune.nue_pid);
const Var kPIDFD = SIMPLEVAR(dune.mvaresult);
const Var kQ = SIMPLEVAR(dune.reco_q);

const HistAxis axis("Reconstructed energy (GeV)",
                    Binning::Simple(40, 0, 10),
                    kRecoE);

const HistAxis axisPID("MVA result",
                       Binning::Simple(60, -1.1, +1.1),
                       kPIDFD);

// One year
const double potND = 1.47e21;
// POT/yr * 3.5yrs * mass correction
const double potFD = 3.5 * 1.47e21 * 40/1.13;

const double pi = 3.14159265359;

const char* stateFname = "nd_fit_state.root";

void Legend(const std::string& title)
{
  TH1* dummy = new TH1F("", "", 1, 0, 1);

  TLegend* leg = new TLegend(.65, .65, .85, .85);
  leg->SetFillStyle(0);
  leg->AddEntry((TObject*)0, ("#bf{"+title+"}").c_str(), "");
  dummy->SetLineColor(kRed);
  leg->AddEntry(dummy->Clone(), "Fake data", "lep");
  dummy->SetLineColor(kBlack);
  leg->AddEntry(dummy->Clone(), "Nominal MC", "l");
  dummy->SetLineColor(kBlue);
  leg->AddEntry(dummy->Clone(), "Best fit", "l");
  dummy->SetLineStyle(kDashed);
  leg->AddEntry(dummy->Clone(), "Best fit #delta = 0", "l");
  leg->Draw();
}

void joint_fit_ND(bool reload = false)
{
  rootlogon(); // style
  Loaders dummyLoaders; // PredictionGenerator insists on this

  osc::IOscCalculatorAdjustable* oscNH = DefaultOscCalc(); // NH, dCP == 0
  osc::IOscCalculatorAdjustable* oscIH = DefaultOscCalcIH(); // IH, dCP == 0

  // all the systematics
  const std::vector<const ISyst*> xsecSysts = GetDUNEXSecSysts(); // uses correlation matrix
  const std::vector<const ISyst*> fluxSysts = GetDUNEFluxSysts(10); // uncorrelated
  const std::vector<const ISyst*> ndSysts = {&kNDEvSyst, &kNDPIDSyst};
  std::vector<const ISyst*> allSysts;
  allSysts.insert( allSysts.end(), xsecSysts.begin(), xsecSysts.end() );
  allSysts.insert( allSysts.end(), fluxSysts.begin(), fluxSysts.end() );
  allSysts.insert( allSysts.end(), ndSysts.begin(), ndSysts.end() );

  std::cout << "Including " << allSysts.size() << " systematics" << std::endl;

  if(reload || TFile(stateFname).IsZombie()){
    SpectrumLoader loaderNDFHC("/dune/data/users/marshalc/CAF_FHC.root");
    SpectrumLoader loaderNDRHC("/dune/data/users/marshalc/CAF_RHC.root");

    auto* loaderNDFHCPOT = loaderNDFHC.LoaderForRunPOT(1);
    auto* loaderNDRHCPOT = loaderNDRHC.LoaderForRunPOT(1);

    SpectrumLoader loaderFDNumuFHC("/pnfs/dune/persistent/TaskForce_AnaTree/far/train/v2.1/numutest.root");
    SpectrumLoader loaderFDNueFHC("/pnfs/dune/persistent/TaskForce_AnaTree/far/train/v2.1/nuetest.root");

    auto* loaderFDNumuFHCBeam  = loaderFDNumuFHC.LoaderForRunPOT(20000001);
    auto* loaderFDNumuFHCNue   = loaderFDNumuFHC.LoaderForRunPOT(20000002);
    auto* loaderFDNumuFHCNuTau = loaderFDNumuFHC.LoaderForRunPOT(20000003);
    auto* loaderFDNumuFHCNC    = loaderFDNumuFHC.LoaderForRunPOT(0);

    auto* loaderFDNueFHCBeam  = loaderFDNueFHC.LoaderForRunPOT(20000001);
    auto* loaderFDNueFHCNue   = loaderFDNueFHC.LoaderForRunPOT(20000002);
    auto* loaderFDNueFHCNuTau = loaderFDNueFHC.LoaderForRunPOT(20000003);
    auto* loaderFDNueFHCNC    = loaderFDNueFHC.LoaderForRunPOT(0);

    SpectrumLoader loaderFDNumuRHC("/pnfs/dune/persistent/TaskForce_AnaTree/far/train/v2.1/anumutest.root");
    SpectrumLoader loaderFDNueRHC("/pnfs/dune/persistent/TaskForce_AnaTree/far/train/v2.1/anuetest.root");

    auto* loaderFDNumuRHCBeam  = loaderFDNumuRHC.LoaderForRunPOT(20000004);
    auto* loaderFDNumuRHCNue   = loaderFDNumuRHC.LoaderForRunPOT(20000005);
    auto* loaderFDNumuRHCNuTau = loaderFDNumuRHC.LoaderForRunPOT(20000006);
    auto* loaderFDNumuRHCNC    = loaderFDNumuRHC.LoaderForRunPOT(0);
    
    auto* loaderFDNueRHCBeam  = loaderFDNueRHC.LoaderForRunPOT(20000004);
    auto* loaderFDNueRHCNue   = loaderFDNueRHC.LoaderForRunPOT(20000005);
    auto* loaderFDNueRHCNuTau = loaderFDNueRHC.LoaderForRunPOT(20000006);
    auto* loaderFDNueRHCNC    = loaderFDNueRHC.LoaderForRunPOT(0);

    NoOscPredictionGenerator genNDFHC(*loaderNDFHCPOT,
                                      axis,
                                      kPIDmu > 0.5 && kQ < 0.);
    PredictionInterp predNDFHC(allSysts,
                               0,
                               genNDFHC,
                               dummyLoaders);

    NoOscPredictionGenerator genNDRHC(*loaderNDRHCPOT,
                                      axis,
                                      kPIDmu > 0.5 && kQ > 0.);
    PredictionInterp predNDRHC(allSysts,
                               0,
                               genNDRHC,
                               dummyLoaders);

    // 0.8 is a random guess at a cut position - should be studied
    DUNENoExtrapPredictionGenerator genFDNumuFHC(*loaderFDNumuFHCBeam, 
                                                 *loaderFDNumuFHCNue,
                                                 *loaderFDNumuFHCNuTau,
                                                 *loaderFDNumuFHCNC,
                                                 axis,
                                                 kPIDFD > 0.8);
    PredictionInterp predFDNumuFHC(allSysts,
                                   oscNH,
                                   genFDNumuFHC,
                                   dummyLoaders);

    // 0.95 is just as arbitrary
    DUNENoExtrapPredictionGenerator genFDNueFHC(*loaderFDNueFHCBeam, 
                                                 *loaderFDNueFHCNue,
                                                 *loaderFDNueFHCNuTau,
                                                 *loaderFDNueFHCNC,
                                                 axis,
                                                 kPIDFD > 0.95);
    PredictionInterp predFDNueFHC(allSysts,
                                  oscNH,
                                  genFDNueFHC,
                                  dummyLoaders);


    DUNENoExtrapPredictionGenerator genFDNumuRHC(*loaderFDNumuRHCBeam, 
                                                 *loaderFDNumuRHCNue,
                                                 *loaderFDNumuRHCNuTau,
                                                 *loaderFDNumuRHCNC,
                                                 axis,
                                                 kPIDFD > 0.8);
    PredictionInterp predFDNumuRHC(allSysts,
                                   oscNH,
                                   genFDNumuRHC,
                                   dummyLoaders);

    DUNENoExtrapPredictionGenerator genFDNueRHC(*loaderFDNueRHCBeam, 
                                                 *loaderFDNueRHCNue,
                                                 *loaderFDNueRHCNuTau,
                                                 *loaderFDNueRHCNC,
                                                 axis,
                                                 kPIDFD > 0.95);
    PredictionInterp predFDNueRHC(allSysts,
                                  oscNH,
                                  genFDNueRHC,
                                  dummyLoaders);


    loaderNDFHC.Go();
    loaderNDRHC.Go();
    loaderFDNumuFHC.Go();
    loaderFDNueFHC.Go();
    loaderFDNumuRHC.Go();
    loaderFDNueRHC.Go();

    TFile fout(stateFname, "RECREATE");
    predNDFHC.SaveTo(fout.mkdir("nd_fhc"));
    predNDRHC.SaveTo(fout.mkdir("nd_rhc"));
    predFDNumuFHC.SaveTo(fout.mkdir("fd_numu_fhc"));
    predFDNueFHC.SaveTo(fout.mkdir("fd_nue_fhc"));
    predFDNumuRHC.SaveTo(fout.mkdir("fd_numu_rhc"));
    predFDNueRHC.SaveTo(fout.mkdir("fd_nue_rhc"));
    std::cout << "Saved state to " << stateFname << std::endl;
  }
  else{
    std::cout << "Loading state from " << stateFname << std::endl;
  }

  TFile fin(stateFname);
  PredictionInterp& predNDFHC = *ana::LoadFrom<PredictionInterp>(fin.GetDirectory("nd_fhc")).release();
  PredictionInterp& predNDRHC = *ana::LoadFrom<PredictionInterp>(fin.GetDirectory("nd_rhc")).release();
  PredictionInterp& predFDNumuFHC = *ana::LoadFrom<PredictionInterp>(fin.GetDirectory("fd_numu_fhc")).release();
  PredictionInterp& predFDNueFHC = *ana::LoadFrom<PredictionInterp>(fin.GetDirectory("fd_nue_fhc")).release();
  PredictionInterp& predFDNumuRHC = *ana::LoadFrom<PredictionInterp>(fin.GetDirectory("fd_numu_rhc")).release();
  PredictionInterp& predFDNueRHC = *ana::LoadFrom<PredictionInterp>(fin.GetDirectory("fd_nue_rhc")).release();
  fin.Close();
  std::cout << "Done loading state" << std::endl;

  // test a shifted fake data set
/*
  DUNEXSecSyst xsshift3(nu_ccqe_3);
  DUNEXSecSyst xsshift2(nu_ccqe_2);
  DUNEXSecSyst xsshift1(nu_ccqe_1);
  std::map<const ISyst*,double> shiftmap;
  shiftmap.emplace(&xsshift3, +3.);
  shiftmap.emplace(&xsshift2, +3.);
  shiftmap.emplace(&xsshift1, +3.);
  shiftmap.emplace(&kNDEvSyst, -1.5);
  SystShifts shifts(shiftmap);

  Spectrum fakeNDFHC = predNDFHC.PredictSyst(0, shifts).FakeData(potND);
  Spectrum fakeNDRHC = predNDRHC.PredictSyst(0, shifts).FakeData(potND);

  Spectrum fakeFDNumuFHC = predFDNumuFHC.PredictSyst(oscNH, shifts).FakeData(potFD);
  Spectrum fakeFDNueFHC = predFDNueFHC.PredictSyst(oscNH, shifts).FakeData(potFD);

  Spectrum fakeFDNumuRHC = predFDNumuRHC.PredictSyst(oscNH, shifts).FakeData(potFD);
  Spectrum fakeFDNueRHC = predFDNueRHC.PredictSyst(oscNH, shifts).FakeData(potFD);
*/

  // Scan 20 values of delta for IH and NH
  // ND predictions will not change with delta CP
  Spectrum fakeNDFHC = predNDFHC.Predict(0).FakeData(potND);
  Spectrum fakeNDRHC = predNDRHC.Predict(0).FakeData(potND);

  const int npoints = 20; // really there will be one more than this
  TGraph * g = new TGraph();
  g->SetName("mcdonals");
  g->SetTitle(";#delta_{CP}/#pi;#Delta#chi^{2}");
  for( int i = 0; i <= npoints; ++i ) {
    double deltaCPpi = -1 + (2*((float)i)/npoints);
    double deltaCP = deltaCPpi*pi; // in units of radians, not pis

    oscNH->SetdCP( deltaCP );
    oscIH->SetdCP( deltaCP );

    Spectrum fakeFDNumuFHC = predFDNumuFHC.Predict(oscNH).FakeData(potFD);
    Spectrum fakeFDNueFHC = predFDNueFHC.Predict(oscNH).FakeData(potFD);
    Spectrum fakeFDNumuRHC = predFDNumuRHC.Predict(oscNH).FakeData(potFD);
    Spectrum fakeFDNueRHC = predFDNueRHC.Predict(oscNH).FakeData(potFD);

    SingleSampleExperiment exptNDFHC(&predNDFHC, fakeNDFHC);
    SingleSampleExperiment exptNDRHC(&predNDRHC, fakeNDRHC);

    SingleSampleExperiment exptFDNumuFHC(&predFDNumuFHC, fakeFDNumuFHC);
    SingleSampleExperiment exptFDNueFHC(&predFDNueFHC, fakeFDNueFHC);

    SingleSampleExperiment exptFDNumuRHC(&predFDNumuRHC, fakeFDNumuRHC);
    SingleSampleExperiment exptFDNueRHC(&predFDNueRHC, fakeFDNueRHC);

    // Joint fit between ND and FD and the covariance matrix
    MultiExperiment expt({&exptNDFHC, &exptNDRHC,
                          &exptFDNumuFHC, &exptFDNueFHC,
                          &exptFDNumuRHC, &exptFDNueRHC,
                          new DUNEXSecCorrelation
                         });

    // Just fit for delta, fix everything else from external data?
    const std::vector<const IFitVar*> oscFitVars = {&kFitDeltaInPiUnits};
    const std::vector<const IFitVar*> emptyOscVars;

    // Set up the fit
    Fitter fit(&expt, oscFitVars, allSysts);
    Fitter fit0(&expt, emptyOscVars, allSysts);
    // Where will we start our search?
    osc::IOscCalculatorAdjustable* oscSeed = DefaultOscCalc();
    SystShifts systSeed = SystShifts::Nominal();
    osc::IOscCalculatorAdjustable* oscSeed0 = DefaultOscCalc();
    SystShifts systSeed0 = SystShifts::Nominal();
    // Do the fit - updates the "seed" variables with the best fit point
    double chi2_floatDelta = fit.Fit(oscSeed, systSeed, Fitter::kQuiet);
    double chi2_delta0 = fit0.Fit(oscSeed0, systSeed0, Fitter::kQuiet);

    std::cout << "True dCP " << deltaCPpi << "pi, best fit " << oscSeed->GetdCP()/pi << "pi, chi2 " << chi2_floatDelta << " delta=0 chi2 " << chi2_delta0 << std::endl;

    g->SetPoint( i, deltaCPpi, chi2_delta0 - chi2_floatDelta );

    // Best fit spectra for all the samples
    Spectrum bfNDFHC = predNDFHC.PredictSyst(0, systSeed);
    Spectrum bfNDRHC = predNDRHC.PredictSyst(0, systSeed);
    Spectrum bfFDNumuFHC = predFDNumuFHC.PredictSyst(oscSeed, systSeed);
    Spectrum bfFDNueFHC = predFDNueFHC.PredictSyst(oscSeed, systSeed);
    Spectrum bfFDNumuRHC = predFDNumuRHC.PredictSyst(oscSeed, systSeed);
    Spectrum bfFDNueRHC = predFDNueRHC.PredictSyst(oscSeed, systSeed);

    Spectrum bfNDFHC0 = predNDFHC.PredictSyst(0, systSeed0);
    Spectrum bfNDRHC0 = predNDRHC.PredictSyst(0, systSeed0);
    Spectrum bfFDNumuFHC0 = predFDNumuFHC.PredictSyst(oscSeed0, systSeed0);
    Spectrum bfFDNueFHC0 = predFDNueFHC.PredictSyst(oscSeed0, systSeed0);
    Spectrum bfFDNumuRHC0 = predFDNumuRHC.PredictSyst(oscSeed0, systSeed0);
    Spectrum bfFDNueRHC0 = predFDNueRHC.PredictSyst(oscSeed0, systSeed0);

    new TCanvas;
    predNDFHC.Predict(0).ToTH1(potND)->Draw("hist");
    bfNDFHC.ToTH1(potND, kBlue)->Draw("hist same");
    fakeNDFHC.ToTH1(potND, kRed)->Draw("ep same");
    Legend("ND FHC");
    gPad->Print( Form("plots/nd_fhc_%1.1fpi.png",deltaCPpi) );

    new TCanvas;
    predNDRHC.Predict(0).ToTH1(potND)->Draw("hist");
    bfNDRHC.ToTH1(potND, kBlue)->Draw("hist same");
    fakeNDRHC.ToTH1(potND, kRed)->Draw("ep same");
    Legend("ND RHC");
    gPad->Print( Form("plots/nd_rhc_%1.1fpi.png",deltaCPpi) );

    new TCanvas;
    predFDNumuFHC.Predict(oscNH).ToTH1(potFD)->Draw("hist");
    bfFDNumuFHC.ToTH1(potFD, kBlue)->Draw("hist same");
    bfFDNumuFHC0.ToTH1(potFD, kBlue, kDashed)->Draw("hist same");
    fakeFDNumuFHC.ToTH1(potFD, kRed)->Draw("ep same");
    Legend("FD FHC #nu_{#mu}");
    gPad->Print( Form("plots/fd_fhc_numu_%1.1fpi.png",deltaCPpi) );

    new TCanvas;
    predFDNueFHC.Predict(oscNH).ToTH1(potFD)->Draw("hist");
    bfFDNueFHC.ToTH1(potFD, kBlue)->Draw("hist same");
    bfFDNueFHC0.ToTH1(potFD, kBlue, kDashed)->Draw("hist same");
    fakeFDNueFHC.ToTH1(potFD, kRed)->Draw("ep same");
    Legend("FD FHC #nu_{e}");
    gPad->Print( Form("plots/fd_fhc_nue_%1.1fpi.png",deltaCPpi) );

    new TCanvas;
    predFDNumuRHC.Predict(oscNH).ToTH1(potFD)->Draw("hist");
    bfFDNumuRHC.ToTH1(potFD, kBlue)->Draw("hist same");
    bfFDNumuRHC0.ToTH1(potFD, kBlue, kDashed)->Draw("hist same");
    fakeFDNumuRHC.ToTH1(potFD, kRed)->Draw("ep same");
    Legend("FD RHC #nu_{#mu}");
    gPad->Print( Form("plots/fd_rhc_numu_%1.1fpi.png",deltaCPpi) );

    new TCanvas;
    predFDNueRHC.Predict(oscNH).ToTH1(potFD)->Draw("hist");
    bfFDNueRHC.ToTH1(potFD, kBlue)->Draw("hist same");
    bfFDNueRHC0.ToTH1(potFD, kBlue, kDashed)->Draw("hist same");
    fakeFDNueRHC.ToTH1(potFD, kRed)->Draw("ep same");
    Legend("FD RHC #nu_{e}");
    gPad->Print( Form("plots/fd_rhc_nue_%1.1fpi.png",deltaCPpi) );
  }

  new TCanvas;
  g->SetMarkerSize(0.5);
  g->Draw("APL");
  gPad->Print( "plots/mcdonalds.png" );


}
