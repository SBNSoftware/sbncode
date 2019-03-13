// Fit including systematics
// cafe demo7.C

#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Vars/FitVarsSterile.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Experiment/SingleSampleExperiment.h"
#include "CAFAna/Analysis/Calcs.h"
#include "CAFAna/Analysis/Surface.h"
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

void demo7()
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

  // Make a vector with all (here only two) the systematics
  std::vector<const ISyst*> allSysts;
  allSysts.push_back(&eSyst);
  allSysts.push_back(&nSyst);

  // List all of the systematics we'll be using
  for(const ISyst* s: allSysts) std::cout << s->ShortName() << "\t\t" << s->LatexName() << std::endl;

  NoExtrapGenerator gen(axEnergy, kSelectionCut);

  PredictionInterp predInterp(allSysts, calc, gen, loaders);

  loaders.Go();

  TCanvas* c1 = new TCanvas("c1");
  TH1* hnom = predInterp.Predict(calc).ToTH1(pot);
  hnom->SetLineWidth(3);
  hnom->Draw("hist");
  hnom->GetYaxis()->SetRangeUser(0, hnom->GetMaximum()*1.3);

  for(int i = 0; i < 100; ++i){
    SystShifts shifts;
    for(const ISyst* s: allSysts) shifts.SetShift(s, gRandom->Gaus());
    predInterp.PredictSyst(calc, shifts).ToTH1(pot, kGray)->Draw("hist same");
  }

  hnom->Draw("hist same");

  c1->SaveAs("demo7_plot1.pdf");

  /// Now show the effect of including the systematic uncertainty in the fit

  c1->Clear();
  c1->SetLogy();
  c1->SetLeftMargin(0.12);
  c1->SetBottomMargin(0.15);

  const Spectrum data = predInterp.Predict(calc).FakeData(pot);
  SingleSampleExperiment expt(&predInterp, data);

  // The regular stats-only contour
  Surface surf(&expt, calc,
               &kFitSinSq2Theta24Sterile, 50, 0, 1,
               &kFitDmSq41Sterile, 75, 0.5, 3);

  surf.DrawBestFit(kBlue);

  // surf.DrawContour(Gaussian68Percent2D(surf), 7, kBlue);
  surf.DrawContour(Gaussian2Sigma2D(surf), kSolid, kBlue);

  gPad->Update();

  // If the list of systematics is long, this fit takes substantially longer.
  // Consider reducing the number of bins
  Surface surfSyst(&expt, calc,
                  &kFitSinSq2Theta24Sterile, 50, 0, 1,
                  &kFitDmSq41Sterile, 75, 0.5, 3,
                  {}, // list of oscillation parameters to profile over
                  allSysts); // list of systs to profile over

  // surfSyst.DrawContour(Gaussian68Percent2D(surfSyst), 7, kRed);
  surfSyst.DrawContour(Gaussian2Sigma2D(surfSyst), kSolid, kRed);

  c1->SaveAs("demo7_plot2.pdf");
}
