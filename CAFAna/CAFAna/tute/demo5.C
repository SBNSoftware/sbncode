// Make oscillated predictions with systematics
// cafe demo5.C

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
      if(sr->sbn.truth.neutrino[0].energy > 2) weight *= 1+0.2*sigma;
      else weight *= 1+0.1*sigma;
    }
  };
  const ToyNormSyst nSyst;

void demo5()
{
  // See demo0.C for explanation of these repeated parts
  const std::string fnameBeam = "/sbnd/app/users/bzamoran/sbncode-v07_11_00/output_largesample_nu_ExampleAnalysis_ExampleSelection.root";
  const std::string fnameSwap = "/sbnd/app/users/bzamoran/sbncode-v07_11_00/output_largesample_oscnue_ExampleAnalysis_ExampleSelection.root";

  // Source of events
  SpectrumLoader loaderBeam(fnameBeam);
  SpectrumLoader loaderSwap(fnameSwap);

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

  PredictionNoExtrap predNom(loaderBeam, loaderSwap, kNullLoader,
                          axEnergy, kSelectionCut);

  // Can set multiple systematics at once like this
  SystShifts bothUp;
  bothUp.SetShift(&eSyst, +1);
  bothUp.SetShift(&nSyst, +1);

  SystShifts bothDn;
  bothDn.SetShift(&eSyst, -1);
  bothDn.SetShift(&nSyst, -1);

  // Each of the constituent OscillatableSpectrum objects within these
  // predictions will be systematically altered as specified in the last
  // argument.
  PredictionNoExtrap predUp(loaderBeam, loaderSwap, kNullLoader,
                            axEnergy, kSelectionCut,
                            bothUp);
  PredictionNoExtrap predDn(loaderBeam, loaderSwap, kNullLoader,
                            axEnergy, kSelectionCut,
                            bothDn);

  // Fill all the nominal and shifted spectra within all three predictions
  loaderBeam.Go();
  loaderSwap.Go();

  const Spectrum sOscNom = predNom.Predict(calc);
  const Spectrum sOscUp = predUp.Predict(calc);
  const Spectrum sOscDn = predDn.Predict(calc);

  TCanvas* c1 = new TCanvas("c1");
  sOscUp.ToTH1(pot, kBlue)->Draw("hist");
  sOscNom.ToTH1(pot)->Draw("hist same");
  sOscDn.ToTH1(pot, kRed)->Draw("hist same");
  c1->SaveAs("demo5_plot1.pdf");
}
