#include "CAFAna/Analysis/Fit.h"
#include "CAFAna/Analysis/CalcsNuFit.h"
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/LoadFromFile.h"
#include "CAFAna/Core/Loaders.h"
#include "CAFAna/Core/Progress.h"
#include "CAFAna/Experiment/MultiExperiment.h"
#include "CAFAna/Experiment/SingleSampleExperiment.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Prediction/PredictionInterp.h"
#include "CAFAna/Prediction/PredictionXSecDiag.h"
#include "CAFAna/Systs/DUNEFluxSysts.h"
#include "CAFAna/Systs/DUNEXSecSysts.h"
#include "CAFAna/Systs/Systs.h"
#include "CAFAna/Vars/FitVars.h"

using namespace ana;

#include "Utilities/rootlogon.C"

#include "OscLib/func/IOscCalculator.h"

#include "StandardRecord/StandardRecord.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2.h"

const Var kRecoE = SIMPLEVAR(dune.Ev_reco);
const Var kPIDFD = SIMPLEVAR(dune.mvaresult);

// 250 MeV bins from 0 to 8GeV
const HistAxis axis("Reconstructed energy (GeV)",
                    Binning::Simple(32, 0, 8),
                    kRecoE);

// POT/yr * 5yrs * mass correction
const double potFD = 5 * 1.47e21 * 40/1.13;

const char* stateFname = "fitter_validation_state.root";
const char* outputFname = "fitter_validation_cafana.root";


void table(FILE* f, IPrediction* p, osc::IOscCalculator* calc)
{
  TH1* hnue = p->PredictComponent(calc, Flavors::kNuMuToNuE, Current::kCC, Sign::kNu).ToTH1(potFD);
  TH1* hnuebar = p->PredictComponent(calc, Flavors::kNuMuToNuE, Current::kCC, Sign::kAntiNu).ToTH1(potFD);
  TH1* hnumu = p->PredictComponent(calc, Flavors::kNuMuToNuMu, Current::kCC, Sign::kNu).ToTH1(potFD);
  TH1* hnumubar = p->PredictComponent(calc, Flavors::kNuMuToNuMu, Current::kCC, Sign::kAntiNu).ToTH1(potFD);
  TH1* hnutau = p->PredictComponent(calc, Flavors::kNuMuToNuTau, Current::kCC, Sign::kNu).ToTH1(potFD);
  TH1* hnutaubar = p->PredictComponent(calc, Flavors::kNuMuToNuTau, Current::kCC, Sign::kAntiNu).ToTH1(potFD);
  TH1* hbeamnue = p->PredictComponent(calc, Flavors::kNuEToNuE, Current::kCC, Sign::kNu).ToTH1(potFD);
  TH1* hbeamnuebar = p->PredictComponent(calc, Flavors::kNuEToNuE, Current::kCC, Sign::kAntiNu).ToTH1(potFD);
  TH1* hnc = p->PredictComponent(calc, Flavors::kAll, Current::kNC, Sign::kBoth).ToTH1(potFD);

  fprintf(f, "E_reco\t\tnue\tnuebar\tnumu\tnumubar\tnutau\ttaubar\tbeamnue\tbeambar\tNC\n");

  for(int i = 1; i <= hnc->GetNbinsX(); ++i){
    fprintf(f, "%g < E < %g\t",
           hnc->GetXaxis()->GetBinLowEdge(i),
           hnc->GetXaxis()->GetBinUpEdge(i));

    fprintf(f, "%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\n",
            hnue->GetBinContent(i),
            hnuebar->GetBinContent(i),
            hnumu->GetBinContent(i),
            hnumubar->GetBinContent(i),
            hnutau->GetBinContent(i),
            hnutaubar->GetBinContent(i),
            hbeamnue->GetBinContent(i),
            hbeamnuebar->GetBinContent(i),
            hnc->GetBinContent(i));
  }
}


double Chisq(IExperiment* expt,
             osc::IOscCalculatorAdjustable* calc,
             bool oscErr, int nfluxErr, int nxsecErr)
{
  if(!oscErr && nfluxErr == 0 && nxsecErr == 0) return expt->ChiSq(calc);

  NuFitPenalizer penalty;

  MultiExperiment exptOscErr({expt, &penalty});

  std::vector<const IFitVar*> oscVars =
    {&kFitDmSq21, &kFitTanSqTheta12,
     &kFitDmSq32Scaled, &kFitSinSqTheta23,
     &kFitTheta13};
  if(!oscErr) oscVars.clear();

  std::vector<const ISyst*> systVars;
  for(int i = 0; i < nfluxErr; ++i) systVars.push_back(GetDUNEFluxSyst(i));
  for(int i = 0; i < nxsecErr; ++i) systVars.push_back(GetDUNEXSecDiagSyst(i));

  Fitter fit(&exptOscErr, oscVars, systVars);

  // Caller doesn't expect all the parameters to get messed up by the fit
  osc::IOscCalculatorAdjustable* calcCopy = calc->Copy();
  SystShifts systSeed;
  const double ret = fit.Fit(calcCopy, systSeed, Fitter::kQuiet);
  delete calcCopy;
  return ret;
}

double ChisqAllCombos(IExperiment* expt, bool oscErr, int nfluxErr, int nxsecErr)
{
  osc::IOscCalculatorAdjustable* oscTest = NuFitOscCalc(+1);
  oscTest->SetdCP(0);
  const double chisq0NH = Chisq(expt, oscTest, oscErr, nfluxErr, nxsecErr);
  oscTest->SetdCP(TMath::Pi());
  const double chisq1NH = Chisq(expt, oscTest, oscErr, nfluxErr, nxsecErr);
  const double chisqNH = std::min(chisq0NH, chisq1NH);

  oscTest = NuFitOscCalc(-1);
  oscTest->SetdCP(0);
  const double chisq0IH = Chisq(expt, oscTest, oscErr, nfluxErr, nxsecErr);
  oscTest->SetdCP(TMath::Pi());
  const double chisq1IH = Chisq(expt, oscTest, oscErr, nfluxErr, nxsecErr);
  const double chisqIH = std::min(chisq0IH, chisq1IH);

  const double chisq = std::min(chisqNH, chisqIH);
  return chisq;
}

void fitter_validation(bool fit = false, bool reload = false)
{
  rootlogon(); // style

  if(reload || TFile(stateFname).IsZombie()){
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

    Loaders dummyLoaders; // PredictionGenerator insists on this

    PredictionNoExtrap predFDNumuFHC(*loaderFDNumuFHCBeam, 
                                     *loaderFDNumuFHCNue,
                                     *loaderFDNumuFHCNuTau,
                                     *loaderFDNumuFHCNC,
                                     axis,
                                     kPIDFD > 0.8);

    PredictionNoExtrap predFDNueFHC(*loaderFDNueFHCBeam, 
                                    *loaderFDNueFHCNue,
                                    *loaderFDNueFHCNuTau,
                                    *loaderFDNueFHCNC,
                                    axis,
                                    kPIDFD > 0.8);

    PredictionNoExtrap predFDNumuRHC(*loaderFDNumuRHCBeam, 
                                     *loaderFDNumuRHCNue,
                                     *loaderFDNumuRHCNuTau,
                                     *loaderFDNumuRHCNC,
                                     axis,
                                     kPIDFD > 0.8);

    PredictionNoExtrap predFDNueRHC(*loaderFDNueRHCBeam, 
                                    *loaderFDNueRHCNue,
                                    *loaderFDNueRHCNuTau,
                                    *loaderFDNueRHCNC,
                                    axis,
                                    kPIDFD > 0.8);


    // Numbers are after sections in Dan's document
    const SystShifts shifts2a(&kNCSyst, +1);

    PredictionNoExtrap predFDNumuFHC2a(*loaderFDNumuFHCBeam, 
                                       *loaderFDNumuFHCNue,
                                       *loaderFDNumuFHCNuTau,
                                       *loaderFDNumuFHCNC,
                                       axis,
                                       kPIDFD > 0.8,
                                       shifts2a);

    PredictionNoExtrap predFDNueFHC2a(*loaderFDNueFHCBeam, 
                                      *loaderFDNueFHCNue,
                                      *loaderFDNueFHCNuTau,
                                      *loaderFDNueFHCNC,
                                      axis,
                                      kPIDFD > 0.8,
                                      shifts2a);

    PredictionNoExtrap predFDNumuRHC2a(*loaderFDNumuRHCBeam, 
                                       *loaderFDNumuRHCNue,
                                       *loaderFDNumuRHCNuTau,
                                       *loaderFDNumuRHCNC,
                                       axis,
                                       kPIDFD > 0.8,
                                       shifts2a);

    PredictionNoExtrap predFDNueRHC2a(*loaderFDNueRHCBeam, 
                                      *loaderFDNueRHCNue,
                                      *loaderFDNueRHCNuTau,
                                      *loaderFDNueRHCNC,
                                      axis,
                                      kPIDFD > 0.8,
                                      shifts2a);

    // TODO - this isn't exactly Dan's prescription
    DUNEFluxSystVector fsv = GetDUNEFluxSysts(256);
    SystShifts shifts2b;
    gRandom->SetSeed(42);
    for(const ISyst* s: fsv) shifts2b.SetShift(s, gRandom->Gaus(0, 1));

    PredictionNoExtrap predFDNumuFHC2b(*loaderFDNumuFHCBeam, 
                                       *loaderFDNumuFHCNue,
                                       *loaderFDNumuFHCNuTau,
                                       *loaderFDNumuFHCNC,
                                       axis,
                                       kPIDFD > 0.8,
                                       shifts2b);

    PredictionNoExtrap predFDNueFHC2b(*loaderFDNueFHCBeam, 
                                      *loaderFDNueFHCNue,
                                      *loaderFDNueFHCNuTau,
                                      *loaderFDNueFHCNC,
                                      axis,
                                      kPIDFD > 0.8,
                                      shifts2b);

    PredictionNoExtrap predFDNumuRHC2b(*loaderFDNumuRHCBeam, 
                                       *loaderFDNumuRHCNue,
                                       *loaderFDNumuRHCNuTau,
                                       *loaderFDNumuRHCNC,
                                       axis,
                                       kPIDFD > 0.8,
                                       shifts2b);

    PredictionNoExtrap predFDNueRHC2b(*loaderFDNueRHCBeam, 
                                      *loaderFDNueRHCNue,
                                      *loaderFDNueRHCNuTau,
                                      *loaderFDNueRHCNC,
                                      axis,
                                      kPIDFD > 0.8,
                                      shifts2b);

    // Flux systematics
    osc::IOscCalculatorAdjustable* inputOsc = NuFitOscCalc(+1);
    DUNENoExtrapPredictionGenerator genFDNumuFHC(*loaderFDNumuFHCBeam, 
                                                 *loaderFDNumuFHCNue,
                                                 *loaderFDNumuFHCNuTau,
                                                 *loaderFDNumuFHCNC,
                                                 axis,
                                                 kPIDFD > 0.8);
    PredictionInterp predFDNumuFHCFlux(GetDUNEFluxSysts(10),
                                       inputOsc,
                                       genFDNumuFHC,
                                       dummyLoaders);

    DUNENoExtrapPredictionGenerator genFDNueFHC(*loaderFDNueFHCBeam, 
                                                 *loaderFDNueFHCNue,
                                                 *loaderFDNueFHCNuTau,
                                                 *loaderFDNueFHCNC,
                                                 axis,
                                                 kPIDFD > 0.8);
    PredictionInterp predFDNueFHCFlux(GetDUNEFluxSysts(10),
                                      inputOsc,
                                      genFDNueFHC,
                                      dummyLoaders);


    DUNENoExtrapPredictionGenerator genFDNumuRHC(*loaderFDNumuRHCBeam, 
                                                 *loaderFDNumuRHCNue,
                                                 *loaderFDNumuRHCNuTau,
                                                 *loaderFDNumuRHCNC,
                                                 axis,
                                                 kPIDFD > 0.8);
    PredictionInterp predFDNumuRHCFlux(GetDUNEFluxSysts(10),
                                       inputOsc,
                                       genFDNumuRHC,
                                       dummyLoaders);

    DUNENoExtrapPredictionGenerator genFDNueRHC(*loaderFDNueRHCBeam, 
                                                 *loaderFDNueRHCNue,
                                                 *loaderFDNueRHCNuTau,
                                                 *loaderFDNueRHCNC,
                                                 axis,
                                                 kPIDFD > 0.8);
    PredictionInterp predFDNueRHCFlux(GetDUNEFluxSysts(10),
                                      inputOsc,
                                      genFDNueRHC,
                                      dummyLoaders);


    // XSec systs

    PredictionXSecDiag predFDNumuFHCXSec(*loaderFDNumuFHCBeam, 
                                         *loaderFDNumuFHCNue,
                                         *loaderFDNumuFHCNuTau,
                                         *loaderFDNumuFHCNC,
                                         axis,
                                         kPIDFD > 0.8);

    PredictionXSecDiag predFDNueFHCXSec(*loaderFDNueFHCBeam, 
                                        *loaderFDNueFHCNue,
                                        *loaderFDNueFHCNuTau,
                                        *loaderFDNueFHCNC,
                                        axis,
                                        kPIDFD > 0.8);

    PredictionXSecDiag predFDNumuRHCXSec(*loaderFDNumuRHCBeam, 
                                         *loaderFDNumuRHCNue,
                                         *loaderFDNumuRHCNuTau,
                                         *loaderFDNumuRHCNC,
                                         axis,
                                         kPIDFD > 0.8);

    PredictionXSecDiag predFDNueRHCXSec(*loaderFDNueRHCBeam, 
                                        *loaderFDNueRHCNue,
                                        *loaderFDNueRHCNuTau,
                                        *loaderFDNueRHCNC,
                                        axis,
                                        kPIDFD > 0.8);


    loaderFDNumuFHC.Go();
    loaderFDNueFHC.Go();
    loaderFDNumuRHC.Go();
    loaderFDNueRHC.Go();

    TFile fout(stateFname, "RECREATE");
    predFDNumuFHC.SaveTo(fout.mkdir("fd_numu_fhc"));
    predFDNueFHC.SaveTo(fout.mkdir("fd_nue_fhc"));
    predFDNumuRHC.SaveTo(fout.mkdir("fd_numu_rhc"));
    predFDNueRHC.SaveTo(fout.mkdir("fd_nue_rhc"));

    predFDNumuFHC2a.SaveTo(fout.mkdir("fd_numu_fhc_2a"));
    predFDNueFHC2a.SaveTo(fout.mkdir("fd_nue_fhc_2a"));
    predFDNumuRHC2a.SaveTo(fout.mkdir("fd_numu_rhc_2a"));
    predFDNueRHC2a.SaveTo(fout.mkdir("fd_nue_rhc_2a"));

    predFDNumuFHC2b.SaveTo(fout.mkdir("fd_numu_fhc_2b"));
    predFDNueFHC2b.SaveTo(fout.mkdir("fd_nue_fhc_2b"));
    predFDNumuRHC2b.SaveTo(fout.mkdir("fd_numu_rhc_2b"));
    predFDNueRHC2b.SaveTo(fout.mkdir("fd_nue_rhc_2b"));

    predFDNumuFHCFlux.SaveTo(fout.mkdir("fd_numu_fhc_flux"));
    predFDNueFHCFlux.SaveTo(fout.mkdir("fd_nue_fhc_flux"));
    predFDNumuRHCFlux.SaveTo(fout.mkdir("fd_numu_rhc_flux"));
    predFDNueRHCFlux.SaveTo(fout.mkdir("fd_nue_rhc_flux"));

    predFDNumuFHCXSec.SaveTo(fout.mkdir("fd_numu_fhc_xsec"));
    predFDNueFHCXSec.SaveTo(fout.mkdir("fd_nue_fhc_xsec"));
    predFDNumuRHCXSec.SaveTo(fout.mkdir("fd_numu_rhc_xsec"));
    predFDNueRHCXSec.SaveTo(fout.mkdir("fd_nue_rhc_xsec"));

    std::cout << "Saved state to " << stateFname << std::endl;
  }
  else{
    std::cout << "Loading state from " << stateFname << std::endl;
  }

  TFile fin(stateFname);
  PredictionNoExtrap& predFDNumuFHC = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_numu_fhc")).release();
  PredictionNoExtrap& predFDNueFHC = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_fhc")).release();
  PredictionNoExtrap& predFDNumuRHC = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_numu_rhc")).release();
  PredictionNoExtrap& predFDNueRHC = *ana::LoadFrom<PredictionNoExtrap>(fin.GetDirectory("fd_nue_rhc")).release();

  IPrediction* predFDNumuFHC2a = ana::LoadFrom<IPrediction>(fin.GetDirectory("fd_numu_fhc_2a")).release();
  IPrediction* predFDNueFHC2a = ana::LoadFrom<IPrediction>(fin.GetDirectory("fd_nue_fhc_2a")).release();
  IPrediction* predFDNumuRHC2a = ana::LoadFrom<IPrediction>(fin.GetDirectory("fd_numu_rhc_2a")).release();
  IPrediction* predFDNueRHC2a = ana::LoadFrom<IPrediction>(fin.GetDirectory("fd_nue_rhc_2a")).release();

  IPrediction* predFDNumuFHC2b = ana::LoadFrom<IPrediction>(fin.GetDirectory("fd_numu_fhc_2b")).release();
  IPrediction* predFDNueFHC2b = ana::LoadFrom<IPrediction>(fin.GetDirectory("fd_nue_fhc_2b")).release();
  IPrediction* predFDNumuRHC2b = ana::LoadFrom<IPrediction>(fin.GetDirectory("fd_numu_rhc_2b")).release();
  IPrediction* predFDNueRHC2b = ana::LoadFrom<IPrediction>(fin.GetDirectory("fd_nue_rhc_2b")).release();

  IPrediction* predFDNumuFHCFlux = ana::LoadFrom<IPrediction>(fin.GetDirectory("fd_numu_fhc_flux")).release();
  IPrediction* predFDNueFHCFlux = ana::LoadFrom<IPrediction>(fin.GetDirectory("fd_nue_fhc_flux")).release();
  IPrediction* predFDNumuRHCFlux = ana::LoadFrom<IPrediction>(fin.GetDirectory("fd_numu_rhc_flux")).release();
  IPrediction* predFDNueRHCFlux = ana::LoadFrom<IPrediction>(fin.GetDirectory("fd_nue_rhc_flux")).release();

  // Has to be explicitly PredictionXSecDiag, otherwise you get back a bare
  // PredictionScaleComp that doesn't do the translation from the diagonalized
  // to underlying space.
  PredictionXSecDiag* predFDNumuFHCXSec = ana::LoadFrom<PredictionXSecDiag>(fin.GetDirectory("fd_numu_fhc_xsec")).release();
  PredictionXSecDiag* predFDNueFHCXSec = ana::LoadFrom<PredictionXSecDiag>(fin.GetDirectory("fd_nue_fhc_xsec")).release();
  PredictionXSecDiag* predFDNumuRHCXSec = ana::LoadFrom<PredictionXSecDiag>(fin.GetDirectory("fd_numu_rhc_xsec")).release();
  PredictionXSecDiag* predFDNueRHCXSec = ana::LoadFrom<PredictionXSecDiag>(fin.GetDirectory("fd_nue_rhc_xsec")).release();

  fin.Close();
  std::cout << "Done loading state" << std::endl;

  TFile* fout = new TFile(outputFname, "RECREATE");

  for(int hie = -1; hie <= +1; hie += 2){
    osc::IOscCalculatorAdjustable* inputOsc = NuFitOscCalc(hie);
    const std::string hieStr = (hie > 0) ? "nh" : "ih";
    for(int deltaIdx = 0; deltaIdx < 4; ++deltaIdx){
      inputOsc->SetdCP(deltaIdx/2.*TMath::Pi());
      const std::string dcpStr = TString::Format("%gpi", deltaIdx/2.).Data();

      FILE* f = fopen(("numu_fhc_"+hieStr+"_"+dcpStr+".txt").c_str(), "w");
      table(f, &predFDNumuFHC, inputOsc);
      fclose(f);

      f = fopen(("nue_fhc_"+hieStr+"_"+dcpStr+".txt").c_str(), "w");
      table(f, &predFDNueFHC, inputOsc);
      fclose(f);

      f = fopen(("numu_rhc_"+hieStr+"_"+dcpStr+".txt").c_str(), "w");
      table(f, &predFDNumuRHC, inputOsc);
      fclose(f);

      f = fopen(("nue_rhc_"+hieStr+"_"+dcpStr+".txt").c_str(), "w");
      table(f, &predFDNueRHC, inputOsc);
      fclose(f);
      
      TH1* hnumufhc = predFDNumuFHC.Predict(inputOsc).ToTH1(potFD);
      TH1* hnuefhc = predFDNueFHC.Predict(inputOsc).ToTH1(potFD);
      TH1* hnumurhc = predFDNumuRHC.Predict(inputOsc).ToTH1(potFD);
      TH1* hnuerhc = predFDNueRHC.Predict(inputOsc).ToTH1(potFD);

      hnumufhc->Write(("numu_fhc_"+hieStr+"_"+dcpStr).c_str());
      hnuefhc->Write(("nue_fhc_"+hieStr+"_"+dcpStr).c_str());
      hnumurhc->Write(("numu_rhc_"+hieStr+"_"+dcpStr).c_str());
      hnuerhc->Write(("nue_rhc_"+hieStr+"_"+dcpStr).c_str());
    } // end for delta
  } // end for hie

  osc::IOscCalculatorAdjustable* inputOsc = NuFitOscCalc(+1);
  inputOsc->SetdCP(0);

  TH1* hnumufhc2a = predFDNumuFHC2a->Predict(inputOsc).ToTH1(potFD);
  TH1* hnuefhc2a = predFDNueFHC2a->Predict(inputOsc).ToTH1(potFD);
  TH1* hnumurhc2a = predFDNumuRHC2a->Predict(inputOsc).ToTH1(potFD);
  TH1* hnuerhc2a = predFDNueRHC2a->Predict(inputOsc).ToTH1(potFD);

  hnumufhc2a->Write("numu_fhc_2a");
  hnuefhc2a->Write("nue_fhc_2a");
  hnumurhc2a->Write("numu_rhc_2a");
  hnuerhc2a->Write("nue_rhc_2a");

  TH1* hnumufhc2b = predFDNumuFHC2b->Predict(inputOsc).ToTH1(potFD);
  TH1* hnuefhc2b = predFDNueFHC2b->Predict(inputOsc).ToTH1(potFD);
  TH1* hnumurhc2b = predFDNumuRHC2b->Predict(inputOsc).ToTH1(potFD);
  TH1* hnuerhc2b = predFDNueRHC2b->Predict(inputOsc).ToTH1(potFD);

  hnumufhc2b->Write("numu_fhc_2b");
  hnuefhc2b->Write("nue_fhc_2b");
  hnumurhc2b->Write("numu_rhc_2b");
  hnuerhc2b->Write("nue_rhc_2b");

  osc::IOscCalculatorAdjustable* osc2e = NuFitOscCalcPlusOneSigma(+1);

  TH1* hnumufhc2e = predFDNumuFHC.Predict(osc2e).ToTH1(potFD);
  TH1* hnuefhc2e = predFDNueFHC.Predict(osc2e).ToTH1(potFD);
  TH1* hnumurhc2e = predFDNumuRHC.Predict(osc2e).ToTH1(potFD);
  TH1* hnuerhc2e = predFDNueRHC.Predict(osc2e).ToTH1(potFD);

  hnumufhc2e->Write("numu_fhc_2e");
  hnuefhc2e->Write("nue_fhc_2e");
  hnumurhc2e->Write("numu_rhc_2e");
  hnuerhc2e->Write("nue_rhc_2e");

  if(!fit) return;

  new TCanvas;
  TH2* axes = new TH2F("", ";#delta_{CP} / #pi;#sigma = #sqrt{#Delta#chi^{2}}", 100, 0, 2, 100, 0, 8);
  axes->Draw();
  
  TGraph* gNH = new TGraph;
  TGraph* gIH = new TGraph;

  TGraph* gNHOscErr = new TGraph;
  TGraph* gIHOscErr = new TGraph;

  TGraph* gNHFlux[10];
  TGraph* gIHFlux[10];
  for(int i = 0; i < 10; ++i){
    gNHFlux[i] = new TGraph;
    gIHFlux[i] = new TGraph;
  }

  TGraph* gNHXSec[10];
  TGraph* gIHXSec[10];
  for(int i = 0; i < 10; ++i){
    gNHXSec[i] = new TGraph;
    gIHXSec[i] = new TGraph;
  }

  for(int hie = -1; hie <= +1; hie += 2){
    Progress prog(hie > 0 ? "NH" : "IH");
    // Chisq explodes at precise CP conservation for some reason
    for(double delta = .001; delta < 2.01; delta += .1){
      osc::IOscCalculatorAdjustable* oscFakeData = NuFitOscCalc(hie);
      oscFakeData->SetdCP(delta*TMath::Pi());

      SingleSampleExperiment exptNueFHC(&predFDNueFHC, predFDNueFHC.Predict(oscFakeData).FakeData(potFD));
      SingleSampleExperiment exptNueRHC(&predFDNueRHC, predFDNueRHC.Predict(oscFakeData).FakeData(potFD));
      SingleSampleExperiment exptNumuFHC(&predFDNumuFHC, predFDNumuFHC.Predict(oscFakeData).FakeData(potFD));
      SingleSampleExperiment exptNumuRHC(&predFDNumuRHC, predFDNumuRHC.Predict(oscFakeData).FakeData(potFD));
      MultiExperiment expt({&exptNueFHC, &exptNueRHC, &exptNumuFHC, &exptNumuRHC});

      SingleSampleExperiment exptNueFHCFlux(predFDNueFHCFlux, predFDNueFHCFlux->Predict(oscFakeData).FakeData(potFD));
      SingleSampleExperiment exptNueRHCFlux(predFDNueRHCFlux, predFDNueRHCFlux->Predict(oscFakeData).FakeData(potFD));
      SingleSampleExperiment exptNumuFHCFlux(predFDNumuFHCFlux, predFDNumuFHCFlux->Predict(oscFakeData).FakeData(potFD));
      SingleSampleExperiment exptNumuRHCFlux(predFDNumuRHCFlux, predFDNumuRHCFlux->Predict(oscFakeData).FakeData(potFD));
      MultiExperiment exptFlux({&exptNueFHCFlux, &exptNueRHCFlux, &exptNumuFHCFlux, &exptNumuRHCFlux});

      SingleSampleExperiment exptNueFHCXSec(predFDNueFHCXSec, predFDNueFHCXSec->Predict(oscFakeData).FakeData(potFD));
      SingleSampleExperiment exptNueRHCXSec(predFDNueRHCXSec, predFDNueRHCXSec->Predict(oscFakeData).FakeData(potFD));
      SingleSampleExperiment exptNumuFHCXSec(predFDNumuFHCXSec, predFDNumuFHCXSec->Predict(oscFakeData).FakeData(potFD));
      SingleSampleExperiment exptNumuRHCXSec(predFDNumuRHCXSec, predFDNumuRHCXSec->Predict(oscFakeData).FakeData(potFD));
      MultiExperiment exptXSec({&exptNueFHCXSec, &exptNueRHCXSec, &exptNumuFHCXSec, &exptNumuRHCXSec});

      const double chisq = ChisqAllCombos(&expt, false, 0, 0);
      const double chisqOscErr = ChisqAllCombos(&expt, true, 0, 0);

      double chisqFlux[10];
      for(int i = 0; i < 10; ++i){
        chisqFlux[i] = ChisqAllCombos(&exptFlux, false, i, 0);
      }

      double chisqXSec[10];
      for(int i = 0; i < 10; ++i){
        chisqXSec[i] = ChisqAllCombos(&exptXSec, false, 0, i);
      }

      if(hie > 0){
        gNH->SetPoint(gNH->GetN(), delta, sqrt(chisq));
        gNHOscErr->SetPoint(gNHOscErr->GetN(), delta, sqrt(chisqOscErr));
        for(int i = 0; i < 10; ++i)
          gNHFlux[i]->SetPoint(gNHFlux[i]->GetN(), delta, sqrt(chisqFlux[i]));
        for(int i = 0; i < 10; ++i)
          gNHXSec[i]->SetPoint(gNHXSec[i]->GetN(), delta, sqrt(chisqXSec[i]));
      }
      else{
        gIH->SetPoint(gIH->GetN(), delta, sqrt(chisq));
        gIHOscErr->SetPoint(gIHOscErr->GetN(), delta, sqrt(chisqOscErr));
        for(int i = 0; i < 10; ++i)
          gIHFlux[i]->SetPoint(gIHFlux[i]->GetN(), delta, sqrt(chisqFlux[i]));
        for(int i = 0; i < 10; ++i)
          gIHXSec[i]->SetPoint(gIHXSec[i]->GetN(), delta, sqrt(chisqXSec[i]));
      }

      prog.SetProgress(delta/2);
    } // end for delta
    prog.Done();
  } // end for hie

  gNH->SetLineColor(kRed);
  gIH->SetLineColor(kBlue);
  gNH->SetLineWidth(2);
  gNH->Draw("l same");
  gIH->SetLineWidth(2);
  gIH->Draw("l same");

  gNHOscErr->SetLineColor(kRed);
  gIHOscErr->SetLineColor(kBlue);
  gNHOscErr->SetLineStyle(7);
  gNHOscErr->SetLineWidth(2);
  gNHOscErr->Draw("l same");
  gIHOscErr->SetLineStyle(7);
  gIHOscErr->SetLineWidth(2);
  gIHOscErr->Draw("l same");

  for(int i = 0; i < 10; ++i){
    gNHFlux[i]->SetLineColor(kRed-7);
    gIHFlux[i]->SetLineColor(kBlue-7);
    gNHFlux[i]->SetLineWidth(1);
    //    gNHFlux[i]->Draw("l same");
    gIHFlux[i]->SetLineWidth(1);
    //    gIHFlux[i]->Draw("l same");
  }

  for(int i = 0; i < 10; ++i){
    gNHXSec[i]->SetLineColor(kRed-7);
    gIHXSec[i]->SetLineColor(kBlue-7);
    gNHXSec[i]->SetLineWidth(1);
    gNHXSec[i]->Draw("l same");
    gIHXSec[i]->SetLineWidth(1);
    gIHXSec[i]->Draw("l same");
  }

  gNH->Write("sens_nh");
  gIH->Write("sens_ih");
  gNHOscErr->Write("sens_nh_oscerr");
  gIHOscErr->Write("sens_ih_oscerr");

  for(int i = 0; i < 10; ++i){
    gNHFlux[i]->Write(TString::Format("sens_nh_flux%d", i).Data());
    gIHFlux[i]->Write(TString::Format("sens_ih_flux%d", i).Data());
  }

  for(int i = 0; i < 10; ++i){
    gNHXSec[i]->Write(TString::Format("sens_nh_xsec%d", i).Data());
    gIHXSec[i]->Write(TString::Format("sens_ih_xsec%d", i).Data());
  }

  std::cout << "Wrote " << outputFname << std::endl;
}
