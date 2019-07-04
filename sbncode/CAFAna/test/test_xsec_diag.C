#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Cuts/TruthCuts.h"

#include "CAFAna/Systs/DUNEXSecSysts.h"
#include "CAFAna/Core/Ratio.h"
#include "CAFAna/Core/SystShifts.h"
#include "CAFAna/Experiment/SingleSampleExperiment.h"
#include "CAFAna/Experiment/MultiExperiment.h"
#include "CAFAna/Prediction/PredictionXSecDiag.h"
#include "CAFAna/Analysis/Calcs.h"
#include "CAFAna/Analysis/Fit.h"
#include "CAFAna/Analysis/Plots.h"
#include "CAFAna/Vars/FitVars.h"

using namespace ana;

#include "Utilities/rootlogon.C"

#include "Utilities/func/MathUtil.h"

#include "StandardRecord/StandardRecord.h"

#include "TCanvas.h"
#include "TGraph.h"
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

const char* stateFname = "test_xsec_diag_state.root";

void test_xsec_diag(bool reload = false)
{
  rootlogon(); // style

  //  LLPerBinFracSystErr::SetError(.01);

  if(reload || TFile(stateFname).IsZombie()){
    SpectrumLoader loaderNDFHC("/dune/data/users/marshalc/NDTF_FGT_FHC_withMEC.root");
    SpectrumLoader loaderNDRHC("/dune/data/users/marshalc/NDTF_FGT_RHC_withMEC.root");

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

    PredictionXSecDiag predNDFHC(*loaderNDFHCPOT, 
                                 axis,
                                 kPIDmu > 0.5 && kQ < 0.);

    PredictionXSecDiag predNDRHC(*loaderNDRHCPOT, 
                                 axis,
                                 kPIDmu > 0.5 && kQ > 0.);

    // 0.8 is a random guess at a cut position - should be studied
    PredictionXSecDiag predFDNumuFHC(*loaderFDNumuFHCBeam, 
                                              *loaderFDNumuFHCNue,
                                              *loaderFDNumuFHCNuTau,
                                              *loaderFDNumuFHCNC,
                                              axis,
                                              kPIDFD > 0.8);

    // As is 0.95
    PredictionXSecDiag predFDNueFHC(*loaderFDNueFHCBeam, 
                                    *loaderFDNueFHCNue,
                                    *loaderFDNueFHCNuTau,
                                    *loaderFDNueFHCNC,
                                    axis,
                                    kPIDFD > 0.95);

    PredictionXSecDiag predFDNumuRHC(*loaderFDNumuRHCBeam, 
                                     *loaderFDNumuRHCNue,
                                     *loaderFDNumuRHCNuTau,
                                     *loaderFDNumuRHCNC,
                                     axis,
                                     kPIDFD > 0.8);

    PredictionXSecDiag predFDNueRHC(*loaderFDNueRHCBeam, 
                                    *loaderFDNueRHCNue,
                                    *loaderFDNueRHCNuTau,
                                    *loaderFDNueRHCNC,
                                    axis,
                                    kPIDFD > 0.95);

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
  //  PredictionXSecDiag& predNDFHC = *PredictionXSecDiag::LoadFrom2(fin.GetDirectory("nd_fhc")).release();

  PredictionXSecDiag& predNDFHC = *ana::LoadFrom<PredictionXSecDiag>(fin.GetDirectory("nd_fhc")).release();
  PredictionXSecDiag& predNDRHC = *ana::LoadFrom<PredictionXSecDiag>(fin.GetDirectory("nd_rhc")).release();
  PredictionXSecDiag& predFDNumuFHC = *ana::LoadFrom<PredictionXSecDiag>(fin.GetDirectory("fd_numu_fhc")).release();
  PredictionXSecDiag& predFDNueFHC = *ana::LoadFrom<PredictionXSecDiag>(fin.GetDirectory("fd_nue_fhc")).release();
  PredictionXSecDiag& predFDNumuRHC = *ana::LoadFrom<PredictionXSecDiag>(fin.GetDirectory("fd_numu_rhc")).release();
  PredictionXSecDiag& predFDNueRHC = *ana::LoadFrom<PredictionXSecDiag>(fin.GetDirectory("fd_nue_rhc")).release();
  fin.Close();
  std::cout << "Done loading state" << std::endl;

  // Make matching FD Asimov fake data, plus some oscillations
  osc::IOscCalculatorAdjustable* inputOsc = DefaultOscCalc();
  inputOsc->SetdCP(1.5*TMath::Pi());

  // What systematic parameters will we shift in the fake data?
  DUNEXSecSyst nushift(nu_MEC_dummy);
  DUNEXSecSyst nubarshift(nubar_MEC_dummy);
  std::map<const ISyst*,double> shiftmap;
  shiftmap.emplace(&nushift, +1.);
  shiftmap.emplace(&nubarshift, +1.);
  SystShifts shifts = shiftmap;

  std::cout << "Truth is " << shifts.ShortName() << std::endl;

  // Make some ND Asimov fake data
  Spectrum fakeNDFHC = predNDFHC.PredictSyst(0, shifts).FakeData(potND);
  Spectrum fakeNDRHC = predNDRHC.PredictSyst(0, shifts).FakeData(potND);

  // We make some FD fake data, but don't bother to use it below
  Spectrum fakeFDNumuFHC = predFDNumuFHC.PredictSyst(inputOsc, shifts).FakeData(potFD);
  Spectrum fakeFDNueFHC = predFDNueFHC.PredictSyst(inputOsc, shifts).FakeData(potFD);

  Spectrum fakeFDNumuRHC = predFDNumuRHC.PredictSyst(inputOsc, shifts).FakeData(potFD);
  Spectrum fakeFDNueRHC = predFDNueRHC.PredictSyst(inputOsc, shifts).FakeData(potFD);

  // Just fit the two ND spectra and make sure we can get out what we put in
  SingleSampleExperiment exptNDFHC(&predNDFHC, fakeNDFHC);
  SingleSampleExperiment exptNDRHC(&predNDRHC, fakeNDRHC);

  MultiExperiment me({&exptNDFHC, &exptNDRHC});

  const std::vector<const ISyst*> systs = GetDUNEXSecDiagSysts();

  for(int nsyst = 1; nsyst <= 32; nsyst *= 2){
    std::vector<const ISyst*> truncSysts = systs;
    truncSysts.resize(nsyst);
    std::cout << std::endl;
    std::cout << "Trying with " << nsyst << " diagonalized systs" << std::endl;
    Fitter fit(&me, {}, truncSysts);
    SystShifts result;
    fit.Fit(0, result);
    std::cout << "Translated:" << std::endl;
    std::cout << predNDFHC.Undiagonalize(result).ShortName() << std::endl;
  }
}
