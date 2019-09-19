#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Cuts/TruthCuts.h"

#include "CAFAna/Systs/DUNEXSecSysts.h"
#include "CAFAna/Systs/DUNEFluxSysts.h"
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

const char* stateFname = "joint_fit_flux_state.root";

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
  leg->Draw();
}

void PlotFluxShiftEffects(osc::IOscCalculator* osc,
                          IPrediction& predNDFHC, IPrediction& predNDRHC,
                          IPrediction& predFDNumuFHC, IPrediction& predFDNueFHC,
                          IPrediction& predFDNumuRHC, IPrediction& predFDNueRHC)
{
  new TCanvas;

  TLegend* leg = new TLegend(.7, .7, .85, .85);
  leg->SetFillStyle(0);
  TH1* dummy = new TH1F("", "", 1, 0, 1);
  dummy->SetLineColor(kRed);
  leg->AddEntry(dummy->Clone(), "+1#sigma", "l");
  dummy->SetLineColor(kBlue);
  leg->AddEntry(dummy->Clone(), "-1#sigma", "l");
  leg->Draw();

  for(int i = 0; i < 10; ++i){
    const ISyst* syst = GetDUNEFluxSyst(i);

    for(bool fd: {false, true}){
      for(bool fhc: {false, true}){
        for(bool nue: {false, true}){
          if(!fd && nue) continue;

          // Triply-nested trinaries. The mark of readable code ;)
          const IPrediction& pred = fd ? (fhc ? (nue ? predFDNueFHC : predFDNumuFHC) : (nue ? predFDNueRHC : predFDNumuRHC)) : (fhc ? predNDFHC : predNDRHC);

          const Spectrum nom = pred.Predict(fd ? osc : 0);

          for(double shift: {-1, +1}){
            const SystShifts ss(syst, shift);

            const Spectrum fake = pred.PredictSyst(fd ? osc : 0, ss);

            const Ratio ratio(fake, nom);
            TH1* h = ratio.ToTH1(shift < 0 ? kBlue : kRed);
            h->Draw(shift < 0 ? "hist" : "hist same");
            h->GetYaxis()->SetRangeUser(0.8, 1.2);
          } // end for shift

          leg->Draw();

          if(fd){
            gPad->Print(TString::Format("plots/flux_syst_ratio_fd_%s_%s_%d.pdf",
                                        nue ? "nue" : "numu",
                                        fhc ? "fhc" : "rhc",
                                        i).Data());
          }
          else{
            gPad->Print(TString::Format("plots/flux_syst_ratio_nd_%s_%d.pdf",
                                        fhc ? "fhc" : "rhc",
                                        i).Data());
          }
        } // end for nue
      } // end for fhc
    } // end for fd
  } // end for i
}

void joint_fit_flux(bool reload = false)
{
  rootlogon(); // style

  Loaders dummyLoaders; // PredictionGenerator insists on this

  osc::IOscCalculatorAdjustable* inputOsc = DefaultOscCalc();
  inputOsc->SetdCP(1.5*TMath::Pi());

  if(reload || TFile(stateFname).IsZombie()){
    SpectrumLoader loaderNDFHC("/dune/data/users/marshalc/NDTF_FGT_FHC.root");
    SpectrumLoader loaderNDRHC("/dune/data/users/marshalc/NDTF_FGT_RHC.root");

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
    PredictionInterp predNDFHC(GetDUNEFluxSysts(10),
                               0,
                               genNDFHC,
                               dummyLoaders);

    NoOscPredictionGenerator genNDRHC(*loaderNDRHCPOT,
                                      axis,
                                      kPIDmu > 0.5 && kQ < 0.);
    PredictionInterp predNDRHC(GetDUNEFluxSysts(10),
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
    PredictionInterp predFDNumuFHC(GetDUNEFluxSysts(10),
                                   inputOsc,
                                   genFDNumuFHC,
                                   dummyLoaders);

    // 0.95 is just as arbitrary
    DUNENoExtrapPredictionGenerator genFDNueFHC(*loaderFDNueFHCBeam, 
                                                 *loaderFDNueFHCNue,
                                                 *loaderFDNueFHCNuTau,
                                                 *loaderFDNueFHCNC,
                                                 axis,
                                                 kPIDFD > 0.95);
    PredictionInterp predFDNueFHC(GetDUNEFluxSysts(10),
                                  inputOsc,
                                  genFDNueFHC,
                                  dummyLoaders);


    DUNENoExtrapPredictionGenerator genFDNumuRHC(*loaderFDNumuRHCBeam, 
                                                 *loaderFDNumuRHCNue,
                                                 *loaderFDNumuRHCNuTau,
                                                 *loaderFDNumuRHCNC,
                                                 axis,
                                                 kPIDFD > 0.8);
    PredictionInterp predFDNumuRHC(GetDUNEFluxSysts(10),
                                   inputOsc,
                                   genFDNumuRHC,
                                   dummyLoaders);

    DUNENoExtrapPredictionGenerator genFDNueRHC(*loaderFDNueRHCBeam, 
                                                 *loaderFDNueRHCNue,
                                                 *loaderFDNueRHCNuTau,
                                                 *loaderFDNueRHCNC,
                                                 axis,
                                                 kPIDFD > 0.95);
    PredictionInterp predFDNueRHC(GetDUNEFluxSysts(10),
                                  inputOsc,
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

  // Makes a lot of plots... best in batch mode
  /*
  predNDFHC.DebugPlots(0, "plots/debug_nd_fhc_%s.pdf");
  predNDRHC.DebugPlots(0, "plots/debug_nd_rhc_%s.pdf");

  predFDNumuFHC.DebugPlots(inputOsc, "plots/debug_fd_numu_fhc_%s.pdf");
  predFDNueFHC.DebugPlots(inputOsc, "plots/debug_fd_nue_fhc_%s.pdf");

  predFDNumuRHC.DebugPlots(inputOsc, "plots/debug_fd_numu_rhc_%s.pdf");
  predFDNueRHC.DebugPlots(inputOsc, "plots/debug_fd_nue_rhc_%s.pdf");

  predNDFHC.DebugPlotsColz(0, "plots/debug_colz_nd_fhc_%s.pdf");
  predNDRHC.DebugPlotsColz(0, "plots/debug_colz_nd_rhc_%s.pdf");

  predFDNumuFHC.DebugPlotsColz(inputOsc, "plots/debug_colz_fd_numu_fhc_%s.pdf");
  predFDNueFHC.DebugPlotsColz(inputOsc, "plots/debug_colz_fd_nue_fhc_%s.pdf");

  predFDNumuRHC.DebugPlotsColz(inputOsc, "plots/debug_colz_fd_numu_rhc_%s.pdf");
  predFDNueRHC.DebugPlotsColz(inputOsc, "plots/debug_colz_fd_nue_rhc_%s.pdf");

  PlotFluxShiftEffects(inputOsc,
                       predNDFHC, predNDRHC,
                       predFDNumuFHC, predFDNueFHC,
                       predFDNumuRHC, predFDNueRHC);
  return;
  */

  // What systematic parameters will we shift in the fake data?
  SystShifts shifts(GetDUNEFluxSyst(0), +1);

  // Make some ND Asimov fake data
  Spectrum fakeNDFHC = predNDFHC.PredictSyst(0, shifts).FakeData(potND);
  Spectrum fakeNDRHC = predNDRHC.PredictSyst(0, shifts).FakeData(potND);

  Spectrum fakeFDNumuFHC = predFDNumuFHC.PredictSyst(inputOsc, shifts).FakeData(potFD);
  Spectrum fakeFDNueFHC = predFDNueFHC.PredictSyst(inputOsc, shifts).FakeData(potFD);

  Spectrum fakeFDNumuRHC = predFDNumuRHC.PredictSyst(inputOsc, shifts).FakeData(potFD);
  Spectrum fakeFDNueRHC = predFDNueRHC.PredictSyst(inputOsc, shifts).FakeData(potFD);

  SingleSampleExperiment exptNDFHC(&predNDFHC, fakeNDFHC);
  SingleSampleExperiment exptNDRHC(&predNDRHC, fakeNDRHC);

  SingleSampleExperiment exptFDNumuFHC(&predFDNumuFHC, fakeFDNumuFHC);
  SingleSampleExperiment exptFDNueFHC(&predFDNueFHC, fakeFDNueFHC);

  SingleSampleExperiment exptFDNumuRHC(&predFDNumuRHC, fakeFDNumuRHC);
  SingleSampleExperiment exptFDNueRHC(&predFDNueRHC, fakeFDNueRHC);

  // Joint fit between ND and FD. No need for covariance matrix with orthoganal
  // flux systs.
  MultiExperiment expt({&exptNDFHC, &exptNDRHC,
                        &exptFDNumuFHC, &exptFDNueFHC,
                        &exptFDNumuRHC, &exptFDNueRHC
                       });

  const std::vector<const IFitVar*> oscFitVars = {&kFitSinSqTheta23,
                                                  &kFitDmSq32Scaled,
                                                  &kFitSinSq2Theta13,
                                                  &kFitDeltaInPiUnits};

  // Use everything - slow
  const std::vector<const ISyst*> systFitVars = GetDUNEFluxSysts(10);

  // Set up the fit
  Fitter fit(&expt, oscFitVars, systFitVars);
  // Where will we start our search?
  osc::IOscCalculatorAdjustable* oscSeed = DefaultOscCalc();
  SystShifts systSeed = SystShifts::Nominal();
  // Do the fit - updates the "seed" variables with the best fit point
  fit.Fit(oscSeed, systSeed);

  // Best fit spectra for all the samples
  Spectrum bfNDFHC = predNDFHC.PredictSyst(0, systSeed);
  Spectrum bfNDRHC = predNDRHC.PredictSyst(0, systSeed);
  Spectrum bfFDNumuFHC = predFDNumuFHC.PredictSyst(oscSeed, systSeed);
  Spectrum bfFDNueFHC = predFDNueFHC.PredictSyst(oscSeed, systSeed);
  Spectrum bfFDNumuRHC = predFDNumuRHC.PredictSyst(oscSeed, systSeed);
  Spectrum bfFDNueRHC = predFDNueRHC.PredictSyst(oscSeed, systSeed);

  // Future work - show the numu/nue/NC sub-components of these

  new TCanvas;
  fakeNDFHC.ToTH1(potND, kRed)->Draw("ep");
  predNDFHC.Predict(0).ToTH1(potND)->Draw("hist same");
  bfNDFHC.ToTH1(potND, kBlue)->Draw("hist same");
  Legend("ND FHC");

  new TCanvas;
  fakeNDRHC.ToTH1(potND, kRed)->Draw("ep");
  predNDRHC.Predict(0).ToTH1(potND)->Draw("hist same");
  bfNDRHC.ToTH1(potND, kBlue)->Draw("hist same");
  Legend("ND RHC");

  new TCanvas;
  fakeFDNumuFHC.ToTH1(potFD, kRed)->Draw("ep");
  predFDNumuFHC.Predict(oscSeed).ToTH1(potFD)->Draw("hist same");
  bfFDNumuFHC.ToTH1(potFD, kBlue)->Draw("hist same");
  Legend("FD FHC #nu_{#mu}");

  new TCanvas;
  fakeFDNueFHC.ToTH1(potFD, kRed)->Draw("ep");
  predFDNueFHC.Predict(oscSeed).ToTH1(potFD)->Draw("hist same");
  bfFDNueFHC.ToTH1(potFD, kBlue)->Draw("hist same");
  Legend("FD FHC #nu_{e}");

  new TCanvas;
  fakeFDNumuRHC.ToTH1(potFD, kRed)->Draw("ep");
  predFDNumuRHC.Predict(oscSeed).ToTH1(potFD)->Draw("hist same");
  bfFDNumuRHC.ToTH1(potFD, kBlue)->Draw("hist same");
  Legend("FD RHC #nu_{#mu}");

  new TCanvas;
  fakeFDNueRHC.ToTH1(potFD, kRed)->Draw("ep");
  predFDNueRHC.Predict(oscSeed).ToTH1(potFD)->Draw("hist same");
  bfFDNueRHC.ToTH1(potFD, kBlue)->Draw("hist same");
  Legend("FD RHC #nu_{e}");
}
