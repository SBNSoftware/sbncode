#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Cuts/TruthCuts.h"

#include "CAFAna/Systs/DUNEXSecSysts.h"
#include "CAFAna/Core/Ratio.h"
#include "CAFAna/Core/SystShifts.h"
#include "CAFAna/Experiment/SingleSampleExperiment.h"
#include "CAFAna/Experiment/MultiExperiment.h"
#include "CAFAna/Analysis/Calcs.h"
#include "CAFAna/Analysis/Fit.h"
#include "CAFAna/Analysis/Plots.h"
#include "CAFAna/Vars/FitVars.h"

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

const char* stateFname = "joint_fit_state.root";

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

void StackPlot(const PredictionScaleComp& pred,
               const std::string& title,
               double pot,
               osc::IOscCalculator* calc = 0)
{
  std::map<EVALORCategory, Color_t> cols;
  cols[nu_ccqe_1] = kRed;
  cols[nu_ccqe_2] = kRed+1;
  cols[nu_ccqe_3] = kRed+2;
  cols[nubar_ccqe_1] = kRed-7;
  cols[nubar_ccqe_2] = kRed-9;
  cols[nubar_ccqe_3] = kRed-10;
  cols[nu_MEC_dummy] = kYellow;
  cols[nubar_MEC_dummy] = kYellow+1;
  cols[nu_cc1piz_1] = kBlue+1;
  cols[nu_cc1piz_2] = kBlue+2;
  cols[nu_cc1piz_3] = kBlue+3;
  cols[nu_cc1pic_1] = kGreen;
  cols[nu_cc1pic_2] = kGreen+1;
  cols[nu_cc1pic_3] = kGreen+2;
  cols[nubar_cc1piz_1] = kBlue-7;
  cols[nubar_cc1piz_2] = kBlue-9;
  cols[nubar_cc1piz_3] = kBlue-10;
  cols[nubar_cc1pic_1] = kGreen-7;
  cols[nubar_cc1pic_2] = kGreen-9;
  cols[nubar_cc1pic_3] = kGreen-10;
  cols[nu_2pi] = kOrange+7;
  cols[nubar_2pi] = kOrange-3;
  cols[nu_dis_1] = kMagenta;
  cols[nu_dis_2] = kMagenta+1;
  cols[nu_dis_3] = kMagenta+2;
  cols[nubar_dis_1] = kMagenta-7;
  cols[nubar_dis_2] = kMagenta-9;
  cols[nubar_dis_3] = kMagenta-10;
  cols[nu_coh] = kCyan;
  cols[nubar_coh] = kCyan+1;
  cols[nu_nc] = kGray+1;
  cols[nubar_nc] = kGray;

  // Draw the total first, so we can check if there are any uncategorized
  // events.
  pred.Predict(calc).ToTH1(pot)->Draw("hist");

  std::vector<std::pair<Spectrum, Color_t>> cats;

  for(int i = 0; i < 32; ++i){
    const EVALORCategory e = EVALORCategory(i);
    Spectrum s = pred.PredictCategory(calc, GetDUNEXSecSyst(e));
    cats.emplace_back(s, cols[e]);
  }
  THStack* stack = ToTHStack(cats, pot);
  stack->Draw("same");
  stack->GetXaxis()->SetTitle("Reconstructed energy (GeV)");
  stack->GetYaxis()->SetTitle("Events");

  TLegend* leg = new TLegend(.7, .3, .9, .85);
  leg->SetFillStyle(0);
  leg->AddEntry((TObject*)0, ("#bf{"+title+"}").c_str(), "");
  for(int i = 0; i < 32; ++i){
    const EVALORCategory e = EVALORCategory(i);
    TH1* dummy = new TH1F("", "", 1, 0, 1);
    dummy->SetFillColor(cols[e]);
    leg->AddEntry(dummy, VALORCategoryName(e).c_str(), "bf");
  }
  leg->Draw();
}

void PlotXSecShiftEffects(osc::IOscCalculator* osc,
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

  for(int i = 0; i < 32; ++i){
    const ISyst* syst = GetDUNEXSecSyst(EVALORCategory(i));

    for(bool fd: {false, true}){
      for(bool fhc: {false, true}){
        for(bool nue: {false, true}){
          if(!fhc && nue) continue;

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
            gPad->Print(TString::Format("plots/syst_ratio_fd_%s_%s_%d.pdf",
                                        nue ? "nue" : "numu",
                                        fhc ? "fhc" : "rhc",
                                        i).Data());
          }
          else{
            gPad->Print(TString::Format("plots/syst_ratio_nd_%s_%d.pdf",
                                        fhc ? "fhc" : "rhc",
                                        i).Data());
          }
        } // end for nue
      } // end for fhc
    } // end for fd
  } // end for i
}

void joint_fit(bool reload = false)
{
  rootlogon(); // style

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

    PredictionScaleComp predNDFHC(*loaderNDFHCPOT, 
                                  axis,
                                  kPIDmu > 0.5 && kQ < 0.,
                                  GetDUNEXSecSysts());

    PredictionScaleComp predNDRHC(*loaderNDRHCPOT, 
                                  axis,
                                  kPIDmu > 0.5 && kQ > 0.,
                                  GetDUNEXSecSysts());

    // 0.8 is a random guess at a cut position - should be studied
    PredictionScaleComp predFDNumuFHC(*loaderFDNumuFHCBeam, 
                                      *loaderFDNumuFHCNue,
                                      *loaderFDNumuFHCNuTau,
                                      *loaderFDNumuFHCNC,
                                      axis,
                                      kPIDFD > 0.8,
                                      GetDUNEXSecSysts());

    PredictionScaleComp predFDNueFHC(*loaderFDNueFHCBeam, 
                                     *loaderFDNueFHCNue,
                                     *loaderFDNueFHCNuTau,
                                     *loaderFDNueFHCNC,
                                     axis,
                                     kPIDFD > 0.95,
                                     GetDUNEXSecSysts());

    PredictionScaleComp predFDNumuRHC(*loaderFDNumuRHCBeam, 
                                      *loaderFDNumuRHCNue,
                                      *loaderFDNumuRHCNuTau,
                                      *loaderFDNumuRHCNC,
                                      axis,
                                      kPIDFD > 0.8,
                                      GetDUNEXSecSysts());

    PredictionScaleComp predFDNueRHC(*loaderFDNueRHCBeam, 
                                     *loaderFDNueRHCNue,
                                     *loaderFDNueRHCNuTau,
                                     *loaderFDNueRHCNC,
                                     axis,
                                     kPIDFD > 0.95,
                                     GetDUNEXSecSysts());

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
  PredictionScaleComp& predNDFHC = *ana::LoadFrom<PredictionScaleComp>(fin.GetDirectory("nd_fhc")).release();
  PredictionScaleComp& predNDRHC = *ana::LoadFrom<PredictionScaleComp>(fin.GetDirectory("nd_rhc")).release();
  PredictionScaleComp& predFDNumuFHC = *ana::LoadFrom<PredictionScaleComp>(fin.GetDirectory("fd_numu_fhc")).release();
  PredictionScaleComp& predFDNueFHC = *ana::LoadFrom<PredictionScaleComp>(fin.GetDirectory("fd_nue_fhc")).release();
  PredictionScaleComp& predFDNumuRHC = *ana::LoadFrom<PredictionScaleComp>(fin.GetDirectory("fd_numu_rhc")).release();
  PredictionScaleComp& predFDNueRHC = *ana::LoadFrom<PredictionScaleComp>(fin.GetDirectory("fd_nue_rhc")).release();
  fin.Close();
  std::cout << "Done loading state" << std::endl;

  // Make matching FD Asimov fake data, plus some oscillations
  osc::IOscCalculatorAdjustable* inputOsc = DefaultOscCalc();
  inputOsc->SetdCP(1.5*TMath::Pi());

  new TCanvas;
  StackPlot(predNDFHC, "ND FHC", potND);
  gPad->Print("cats_nd_fhc_numu.pdf");
  new TCanvas;
  StackPlot(predNDRHC, "ND RHC", potND);
  gPad->Print("cats_nd_rhc_numu.pdf");

  new TCanvas;
  StackPlot(predFDNumuFHC, "FD FHC #nu_{#mu}", potFD, inputOsc);
  gPad->Print("cats_fd_fhc_numu.pdf");
  new TCanvas;
  StackPlot(predFDNueFHC,  "FD FHC #nu_{e}",   potFD, inputOsc);
  gPad->Print("cats_fd_fhc_nue.pdf");
  new TCanvas;
  StackPlot(predFDNumuRHC, "FD RHC #nu_{#mu}", potFD, inputOsc);
  gPad->Print("cats_fd_rhc_numu.pdf");
  new TCanvas;
  StackPlot(predFDNueRHC,  "FD RHC #nu_{e}",   potFD, inputOsc);
  gPad->Print("cats_fd_rhc_nue.pdf");

  PlotXSecShiftEffects(inputOsc,
                       predNDFHC, predNDRHC,
                       predFDNumuFHC, predFDNueFHC,
                       predFDNumuRHC, predFDNueRHC);

  // What systematic parameters will we shift in the fake data?
  DUNEXSecSyst nushift(nu_MEC_dummy);
  DUNEXSecSyst nubarshift(nubar_MEC_dummy);
  std::map<const ISyst*,double> shiftmap;
  shiftmap.emplace(&nushift, +1.);
  shiftmap.emplace(&nubarshift, +1.);
  SystShifts shifts(shiftmap);

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

  // Joint fit between ND and FD and the covariance matrix
  MultiExperiment expt({&exptNDFHC, &exptNDRHC,
                        &exptFDNumuFHC, &exptFDNueFHC,
                        &exptFDNumuRHC, &exptFDNueRHC,
                        new DUNEXSecCorrelation
                       });

  const std::vector<const IFitVar*> oscFitVars = {&kFitSinSqTheta23,
                                                  &kFitDmSq32Scaled,
                                                  &kFitSinSq2Theta13,
                                                  &kFitDeltaInPiUnits};

  // Use everything - slow
  const std::vector<const ISyst*> systFitVars = GetDUNEXSecSysts();
  // Fit a reduced list of variables of interest
  //    const std::vector<const ISyst*> systFitVars = {GetDUNEXSecSyst(nu_MEC_dummy),
  //                                                   GetDUNEXSecSyst(nubar_MEC_dummy)};

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
