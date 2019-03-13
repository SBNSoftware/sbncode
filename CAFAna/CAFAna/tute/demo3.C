// Make a simple contour
// cafe demo3.C

#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Analysis/Calcs.h"
#include "OscLib/func/OscCalculatorSterile.h"
#include "StandardRecord/StandardRecord.h"
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

void demo3()
{
  // See demo0.C for explanation of these repeated parts
  const std::string fnameBeam = "/sbnd/app/users/bzamoran/sbncode-v07_11_00/output_largesample_nu_ExampleAnalysis_ExampleSelection.root";
  const std::string fnameSwap = "/sbnd/app/users/bzamoran/sbncode-v07_11_00/output_largesample_oscnue_ExampleAnalysis_ExampleSelection.root";

  // Source of events
  SpectrumLoader loaderBeam(fnameBeam);
  SpectrumLoader loaderSwap(fnameSwap);

  // And now add Icarus data
  const std::string fnameBeam2 = "/sbnd/app/users/bzamoran/sbncode-v07_11_00/output_icarus_nu_ExampleAnalysis_ExampleSelection.root";
  const std::string fnameSwap2 = "/sbnd/app/users/bzamoran/sbncode-v07_11_00/output_icarus_oscnue_ExampleAnalysis_ExampleSelection.root";

  // Source of events
  SpectrumLoader loaderBeam2(fnameBeam2);
  SpectrumLoader loaderSwap2(fnameSwap2);

  const Var kRecoEnergy({}, // ToDo: smear with some resolution
                        [](const caf::StandardRecord* sr)
                        {
                          double fE = sr->sbn.truth.neutrino[0].energy;
                          TRandom3 r(floor(fE*10000));
                          double smear = r.Gaus(1, 0.05); // Flat 5% E resolution
                          return fE*smear;
                        });

  const Binning binsEnergy = Binning::Simple(50, 0, 5);
  const HistAxis axEnergy("Fake reconsturcted energy (GeV)", binsEnergy, kRecoEnergy);

  // Fake POT: we need to sort this out in the files first
  const double pot = 6.e20;

  const Cut kSelectionCut({},
                       [](const caf::StandardRecord* sr)
                       {
                         double fE = sr->sbn.truth.neutrino[0].energy;
                         TRandom3 r(floor(fE*10000));
                         bool isCC = sr->sbn.truth.neutrino[0].iscc;
                         double p = r.Uniform();
                         // 80% eff for CC, 10% for NC
                         if(isCC) return p < 0.8;
                         else return p < 0.10;
                       });

  PredictionNoExtrap pred(loaderBeam, loaderSwap, kNullLoader,
                          axEnergy, kSelectionCut);

  PredictionNoExtrap pred2(loaderBeam2, loaderSwap2, kNullLoader,
                          axEnergy, kSelectionCut);

  loaderBeam.Go();
  loaderSwap.Go();

  loaderBeam2.Go();
  loaderSwap2.Go();

  // Calculator
  osc::OscCalculatorSterile* calc = DefaultSterileCalc(4);
  calc->SetL(0.11); // SBND only, temporary
  calc->SetAngle(2, 4, 0.55);
  calc->SetDm(4, 1); // Some dummy values

  TMarker* trueValues = new TMarker(pow(TMath::Sin(2*calc->GetAngle(2,4)),2), calc->GetDm(4), kFullCircle);
  trueValues->SetMarkerColor(kRed);

  // To make a fit we need to have a "data" spectrum to compare to our MC
  // Prediction object
  const Spectrum data = pred.Predict(calc).FakeData(pot);
  SingleSampleExperiment expt(&pred, data);

  // A Surface evaluates the experiment's chisq across a grid
  Surface surf(&expt, calc,
               &kFitSinSq2Theta24Sterile, 50, 0, 1,
               &kFitDmSq41Sterile, 75, 0.5, 3);

  TCanvas* c1 = new TCanvas("c1");
  c1->SetLogy();
  c1->SetLeftMargin(0.12);
  c1->SetBottomMargin(0.15);
  //surf.Draw();
  surf.DrawBestFit(kBlue);
  trueValues->Draw();

  // In a full Feldman-Cousins analysis you need to provide a critical value
  // surface to be able to draw a contour. But we provide these helper
  // functions to use the gaussian up-values.
  TH2* crit1sig = Gaussian68Percent2D(surf);
  TH2* crit2sig = Gaussian2Sigma2D(surf);

  surf.DrawContour(crit1sig, 7, kBlue);
  surf.DrawContour(crit2sig, kSolid, kBlue);

  c1->SaveAs("demo3_plot1.pdf");

  // Let's now try adding a second experiment, for instance Icarus
  calc->SetL(0.6); // Turn of Icarus
  const Spectrum data2 = pred2.Predict(calc).FakeData(pot);
  SingleSampleExperiment expt2(&pred2, data2);

  MultiExperimentSBN multiExpt({&expt, &expt2}, {0.11, 0.6});

  Surface surf2(&expt2, calc,
               &kFitSinSq2Theta24Sterile, 50, 0, 1,
               &kFitDmSq41Sterile, 75, 0.5, 3);

  c1->Clear(); // just in case
  surf2.DrawBestFit(kGreen+2);
  trueValues->Draw();
  TH2* crit1sig2 = Gaussian68Percent2D(surf2);
  TH2* crit2sig2 = Gaussian2Sigma2D(surf2);

  surf2.DrawContour(crit1sig2, 7, kGreen+2);
  surf2.DrawContour(crit2sig2, kSolid, kGreen+2);

  c1->SaveAs("demo3_plot2.pdf");

  Surface surfMulti(&multiExpt, calc,
               &kFitSinSq2Theta24Sterile, 50, 0, 1,
               &kFitDmSq41Sterile, 75, 0.5, 3);

  c1->Clear(); // just in case
  surfMulti.DrawBestFit(kBlue);
  trueValues->Draw();
  TH2* crit1sigMulti = Gaussian68Percent2D(surfMulti);
  TH2* crit2sigMulti = Gaussian2Sigma2D(surfMulti);

  surfMulti.DrawContour(crit1sigMulti, 7, kBlue);
  surfMulti.DrawContour(crit2sigMulti, kSolid, kBlue);

  c1->SaveAs("demo3_plot3.pdf");

  // // We can now fit including profiling over the (relevant) unknown 
  // // oscillation parameters. We can also include constraints for those.
  // // ToDo: the complete general case, with multiple experiments and profiling

  // // Experimental constraints (NuFit 4.0)
  // GaussianConstraint th23Constraint(&kFitTheta23InDegreesSterile, 49.7, 1.1);
  // GaussianConstraint dmsq32Constraint(&kFitDmSq32Sterile, 2.53e-3, 0.03e-3);
  // MultiExperiment multiExpt3({&th23Constraint, &dmsq32Constraint, &expt});

  // std::vector< const IFitVar * > kProfiledVars = {&kFitTheta23InDegreesSterile, 
  //                                                 &kFitDmSq32Sterile,
  //                                                 &kFitDelta14InPiUnitsSterile,
  //                                                 &kFitDelta24InPiUnitsSterile};

  // Surface surf3(&multiExpt3, calc,
  //              &kFitSinSq2Theta24Sterile, 50, 0, 1,
  //              &kFitDmSq41Sterile, 75, 0.5, 3,
  //              kProfiledVars,{});

  // c1->Clear(); // just in case
  // surf3.DrawBestFit(kBlue);
  // trueValues->Draw();

  // TH2* crit1sig3 = Gaussian68Percent2D(surf3);
  // TH2* crit2sig3 = Gaussian2Sigma2D(surf3);

  // surf3.DrawContour(crit1sig3, 7, kBlue);
  // surf3.DrawContour(crit2sig3, kSolid, kBlue);

  // c1->SaveAs("demo3_plot4.pdf");

}
