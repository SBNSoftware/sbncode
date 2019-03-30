#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Analysis/Calcs.h"
#include "OscLib/func/OscCalculatorSterile.h"
#include "StandardRecord/Proxy/SRProxy.h"
#include "TCanvas.h"
#include "TH1.h"
#include "CAFAna/Vars/FitVarsSterile.h"

// New includes
#include "CAFAna/Analysis/Surface.h"
#include "CAFAna/Experiment/SingleSampleExperiment.h"
#include "CAFAna/Experiment/MultiExperiment.h"
#include "CAFAna/Experiment/MultiExperimentSBN.h"
#include "CAFAna/Experiment/GaussianConstraint.h"

// Random numbers to fake an efficiency and resolution
#include "TRandom3.h"

#include "TMarker.h"

using namespace ana;

void numu()
{
  // See demo0.C for explanation of these repeated parts
  const std::string fnameBeam = "/pnfs/sbnd/persistent/users/gputnam/numu_simulation_12_05_2018/processed/output_SBNOsc_NumuSelection_Modern_SBND.root";
  // const std::string fnameSwap = "/sbnd/app/users/bzamoran/sbncode-v07_11_00/output_largesample_oscnue_ExampleAnalysis_ExampleSelection.root";

  // Source of events
  SpectrumLoader loaderBeam(fnameBeam);
  // SpectrumLoader loaderSwap(fnameSwap);

  // And now add Icarus data
  const std::string fnameBeam2 = "/pnfs/sbnd/persistent/users/gputnam/numu_simulation_12_05_2018/processed/output_SBNOsc_NumuSelection_Modern_Icarus.root";
  // const std::string fnameSwap2 = "/sbnd/app/users/bzamoran/sbncode-v07_11_00/output_icarus_oscnue_ExampleAnalysis_ExampleSelection.root";

  // Source of events
  SpectrumLoader loaderBeam2(fnameBeam2);
  // SpectrumLoader loaderSwap2(fnameSwap2);

  const double sbndPOT = 6.6e20 * 1.e20 / 47883366000000000000.000000;
  const double icarusPOT = 6.6e20 * 1.e20 / 118884300000000000000.000000;

  const Var kRecoEnergy(// ToDo: smear with some resolution
                        [](const caf::SRProxy* sr)
                        {
                          double fE = sr->sbn.truth.neutrino[0].energy;
                          TRandom3 r(floor(fE*10000));
                          double smear = r.Gaus(1, 0.05); // Flat 5% E resolution
                          return fE*smear;
                        });

  const Binning binsEnergy = Binning::Simple(30, 0, 3);
  const HistAxis axEnergy("Fake reconsturcted energy (GeV)", binsEnergy, kRecoEnergy);

  const Cut kSelectionCut([](const caf::SRProxy* sr)
                          {
                            double fE = sr->sbn.truth.neutrino[0].energy;
                            TRandom3 r(floor(fE*10000));
                            bool isCC = sr->sbn.truth.neutrino[0].iscc;
                            double p = r.Uniform();
                            // 80% eff for CC, 10% for NC
                            if(isCC) return p < 0.8;
                            else return p < 0.10;
                          });

  PredictionNoExtrap pred(loaderBeam, kNullLoader, kNullLoader,
                          axEnergy, kSelectionCut);

  PredictionNoExtrap pred2(loaderBeam2, kNullLoader, kNullLoader,
                          axEnergy, kSelectionCut);

  loaderBeam.Go();
  // loaderSwap.Go();

  loaderBeam2.Go();
  // loaderSwap2.Go();

  // Calculator
  osc::OscCalculatorSterile* calc = DefaultSterileCalc(4);
  calc->SetAngle(2, 4, 0);
  calc->SetDm(4, 0); // Make these explicit
  calc->SetL(0.11); // SBND only, temporary
  calc->SetAngle(2, 3, M_PI/4);
  calc->SetDm(3, 2.45e-3); // make these explicit

  // TMarker* trueValues = new TMarker(pow(TMath::Sin(2*calc->GetAngle(2,4)),2), calc->GetDm(4), kFullCircle);
  // trueValues->SetMarkerColor(kRed);

  // To make a fit we need to have a "data" spectrum to compare to our MC
  // Prediction object
  const Spectrum data = pred.Predict(calc).FakeData(sbndPOT);
  SingleSampleExperiment expt(&pred, data);

  TFile* fOutput = new TFile("4flavour/Surfaces.root","RECREATE");

  // A Surface evaluates the experiment's chisq across a grid
  Surface surf(&expt, calc,
               &kFitSinSqTheta23Sterile, 50, 0, 1,
               &kFitDmSq32Sterile, 50, 0, 5e-3);

  surf.SaveTo(fOutput->mkdir("surf"));

  TCanvas* c1 = new TCanvas("c1");
  // c1->SetLogx();
  // c1->SetLogy();
  c1->SetLeftMargin(0.12);
  c1->SetBottomMargin(0.15);
  //surf.Draw();
  // surf.DrawBestFit(kBlue);
  // trueValues->Draw();

  // In a full Feldman-Cousins analysis you need to provide a critical value
  // surface to be able to draw a contour. But we provide these helper
  // functions to use the gaussian up-values.
  TH2* crit1sig = Gaussian68Percent2D(surf);
  TH2* crit2sig = Gaussian2Sigma2D(surf);

  surf.DrawContour(crit1sig, 7, kBlue);
  surf.DrawContour(crit2sig, kSolid, kBlue);

  c1->SaveAs("4flavour/numu_plot1.pdf");

  // Let's now try adding a second experiment, for instance Icarus
  calc->SetL(0.6); // Turn of Icarus
  const Spectrum data2 = pred2.Predict(calc).FakeData(icarusPOT);
  SingleSampleExperiment expt2(&pred2, data2);

  MultiExperimentSBN multiExpt({&expt, &expt2}, {0.11, 0.6});

  Surface surf2(&expt2, calc,
               &kFitSinSqTheta23Sterile, 50, 0, 1,
               &kFitDmSq32Sterile, 50, 0, 5e-3);

  surf2.SaveTo(fOutput->mkdir("surf2"));

  c1->Clear(); // just in case
  // surf2.DrawBestFit(kGreen+2);
  // trueValues->Draw();
  TH2* crit1sig2 = Gaussian68Percent2D(surf2);
  TH2* crit2sig2 = Gaussian2Sigma2D(surf2);

  surf2.DrawContour(crit1sig2, 7, kGreen+2);
  surf2.DrawContour(crit2sig2, kSolid, kGreen+2);

  c1->SaveAs("4flavour/numu_plot2.pdf");

  Surface surfMulti(&multiExpt, calc,
               &kFitSinSqTheta23Sterile, 50, 0, 1,
               &kFitDmSq32Sterile, 50, 0, 5e-3);

  surfMulti.SaveTo(fOutput->mkdir("surfMulti"));

  c1->Clear(); // just in case
  // surfMulti.DrawBestFit(kBlue);
  // trueValues->Draw();
  TH2* crit1sigMulti = Gaussian68Percent2D(surfMulti);
  TH2* crit2sigMulti = Gaussian2Sigma2D(surfMulti);

  surfMulti.DrawContour(crit1sigMulti, 7, kBlue);
  surfMulti.DrawContour(crit2sigMulti, kSolid, kBlue);

  c1->SaveAs("4flavour/numu_plot3.pdf");

  fOutput->Close();

}
