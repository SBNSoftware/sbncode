// Introduces PredictionInterp
// cafe demo6.C

#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Analysis/Calcs.h"
#include "StandardRecord/StandardRecord.h"
#include "OscLib/func/OscCalculatorSterile.h"
#include "TCanvas.h"
#include "TH1.h"

// new includes
#include "CAFAna/Core/Loaders.h"
#include "CAFAna/Prediction/PredictionInterp.h"
#include "CAFAna/Prediction/PredictionGenerator.h"

// Random numbers to fake an efficiency and resolution
#include "TRandom3.h"

using namespace ana;

class ToyEnergyScaleSyst: public ISyst
  {
  public:
    ToyEnergyScaleSyst() : ISyst("toyEScale", "Toy Energy Scale") {}
    void Shift(double sigma,
               Restorer& restore,
               caf::StandardRecord* sr,
               double& weight) const override
    {
      restore.Add(sr->sbn.truth.neutrino[0].energy);
      const double scale = 1 + .03*sigma; // 3% E scale syst.
      sr->sbn.truth.neutrino[0].energy *= scale;
    }
  };
  const ToyEnergyScaleSyst eSyst;

class ToyNormSyst: public ISyst
  {
  public:
    ToyNormSyst() : ISyst("toyNorm", "Toy Norm Scale") {}
    void Shift(double sigma,
               Restorer& restore,
               caf::StandardRecord* sr,
               double& weight) const override
    {
      if(sr->sbn.truth.neutrino[0].energy > 2) weight *= TMath::Max(0., 1+0.2*sigma);
      else weight *= TMath::Max(0., 1+0.1*sigma);
    }
  };
  const ToyNormSyst nSyst;

void demo6()
{
  // See demo0.C for explanation of these repeated parts
  const std::string fnameBeam = "/sbnd/app/users/bzamoran/sbncode-v07_11_00/output_largesample_nu_ExampleAnalysis_ExampleSelection.root";
  const std::string fnameSwap = "/sbnd/app/users/bzamoran/sbncode-v07_11_00/output_largesample_oscnue_ExampleAnalysis_ExampleSelection.root";

  // Source of events
  Loaders loaders;

  loaders.SetLoaderPath( fnameBeam,  caf::kFARDET,  Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);
  loaders.SetLoaderPath( fnameSwap,  caf::kFARDET,  Loaders::kMC,   ana::kBeam, Loaders::kFluxSwap);

  const Var kRecoEnergy({}, // ToDo: smear with some resolution
                        [](const caf::StandardRecord* sr)
                        {
                          double fE = sr->sbn.truth.neutrino[0].energy;
                          TRandom3 r(floor(fE*10000));
                          double smear = r.Gaus(1, 0.05); // Flat 5% E resolution
                          return fE*smear;
                        });

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

  const Binning binsEnergy = Binning::Simple(50, 0, 5);
  const HistAxis axEnergy("Fake reconsturcted energy (GeV)", binsEnergy, kRecoEnergy);

  // Fake POT: we need to sort this out in the files first
  const double pot = 6.e20;

  // Calculator
  osc::OscCalculatorSterile* calc = DefaultSterileCalc(4);
  calc->SetL(0.11); // SBND only, temporary
  calc->SetAngle(2, 4, 0.55);
  calc->SetDm(4, 1); // Some dummy values

  // We're going to use a PredictionInterp that will allow us to interpolate to
  // any values of the systematic parameters. Internally that works by creating
  // various predictions at different values of the paramters, so we need to
  // add this extra layer of indirection to allow it to create those.
  NoExtrapGenerator gen(axEnergy, kSelectionCut);

  // PredictionInterp needs:
  // - The list of systematics it should be ready to interpolate over
  // - A "seed" oscillation calculator it will do its expansions around
  // - The generator from above
  // - The dummy loaders (sorry)

  PredictionInterp predInterp({&eSyst, &nSyst},
                              calc,
                              gen,
                              loaders);

  // Fill all the different variants of the predictions that PredictionInterp
  // needs to make.
  loaders.Go();

  // Make some nice plots of what the interpolation looks like in each bin
  predInterp.DebugPlots(calc);

  // We can generate predictions at whatever values of the systematic shifts we
  // want. Prove it by plotting 100 different possible "universes".
  TCanvas* c1 = new TCanvas("c1");

  TH1* hnom = predInterp.Predict(calc).ToTH1(pot);
  hnom->SetLineWidth(3);
  hnom->Draw("hist");
  hnom->GetYaxis()->SetRangeUser(0, hnom->GetMaximum()*1.3);

  for(int i = 0; i < 100; ++i){
    SystShifts s;
    s.SetShift(&eSyst, gRandom->Gaus());
    s.SetShift(&nSyst, gRandom->Gaus());

    TH1* hshift = predInterp.PredictSyst(calc, s).ToTH1(pot);
    hshift->SetLineColor(kGray);
    hshift->Draw("hist same");
  }

  hnom->Draw("hist same");
  c1->SaveAs("demo6_plot1.pdf");
}
