// Make oscillated predictions
// cafe demo1.C

#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
// #include "CAFAna/Cuts/TruthCuts.h"
#include "StandardRecord/StandardRecord.h"
#include "TCanvas.h"
#include "TH1.h"

// New includes for this macro
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Analysis/Calcs.h"
#include "OscLib/func/OscCalculatorSterile.h"

// Random numbers to fake an efficiency and resolution
#include "TRandom3.h"

using namespace ana;

void demo1()
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

  // In many cases it's easier to form them from existing Vars like this
  //  const Cut kPassesMVA = kMVANumu > 0;

  // A Prediction is an objects holding a variety of "OscillatableSpectrum"
  // objects, one for each original and final flavour combination.
  PredictionNoExtrap pred(loaderBeam, loaderSwap, kNullLoader,
                          axEnergy, kSelectionCut);

  // This call will fill all of the constituent parts of the prediction
  loaderBeam.Go();
  loaderSwap.Go();

  // We can extract a total prediction unoscillated
  const Spectrum sUnosc = pred.PredictUnoscillated();
  // Or oscillated, in this case using reasonable parameters from
  // Analysis/Calcs.h
  osc::OscCalculatorSterile* calc = DefaultSterileCalc(4);
  calc->SetL(0.11); // SBND only, temporary
  calc->SetAngle(2, 4, 0.55);
  calc->SetDm(4, 1); // Some dummy values

  const Spectrum sOsc = pred.Predict(calc);

  // And we can break things down by flavour
  const Spectrum sUnoscNC = pred.PredictComponent(calc,
                                                  Flavors::kAll,
                                                  Current::kNC,
                                                  Sign::kBoth);

  // Plot what we have so far
  TCanvas* c1 = new TCanvas("c1");
  sUnosc.ToTH1(pot)->Draw("hist");
  sUnoscNC.ToTH1(pot, kBlue)->Draw("hist same");
  sOsc.ToTH1(pot, kRed)->Draw("hist same");
  c1->SaveAs("demo1_plot1.pdf");

  // "Fake" data is synonymous with the Asimov data sample
  sOsc.ToTH1(pot, kRed)->Draw("hist");
  sUnoscNC.ToTH1(pot, kBlue)->Draw("hist same");
  sOsc.FakeData(pot).ToTH1(pot)->Draw("ep same");
  c1->SaveAs("demo1_plot2.pdf");

  // While "mock" data has statistical fluctuations in
  sOsc.ToTH1(pot, kRed)->Draw("hist");
  sUnoscNC.ToTH1(pot, kBlue)->Draw("hist same");
  sOsc.MockData(pot).ToTH1(pot)->Draw("ep same");
  c1->SaveAs("demo1_plot3.pdf");
}
