// Make a simple spectrum plot with systematic shifts
// cafe demo4.C

#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"

#include "CAFAna/Cuts/TruthCuts.h"

#include "StandardRecord/StandardRecord.h"

#include "TCanvas.h"
#include "TH1.h"
#include "TPad.h"

// Random numbers to fake an efficiency and resolution
#include "TRandom3.h"

using namespace ana;

void demo4()
{
  // See demo0.C for explanation of these repeated parts
  const std::string fnameBeam = "/sbnd/app/users/bzamoran/sbncode-v07_11_00/output_largesample_nu_ExampleAnalysis_ExampleSelection.root";

  // Source of events
  SpectrumLoader loaderBeam(fnameBeam);

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

  // This is the nominal energy spectrum
  Spectrum sEnergy(loaderBeam, axEnergy, kSelectionCut);

  // Systematics work by modifying the event record before it's filled into the
  // spectrum. These generally should be added into a header like Systs.h
  // rather than implemented inline like this.
  class ToyEnergyScaleSyst: public ISyst
  {
  public:
    ToyEnergyScaleSyst() : ISyst("toyEScale", "Toy Energy Scale") {}

    // Function that will be called to actually do the shift
    void Shift(double sigma,
               Restorer& restore,
               caf::StandardRecord* sr,
               double& weight) const override
    {
      // First - register all the variables that will need to be restored to
      // return the record to nominal
      restore.Add(sr->sbn.truth.neutrino[0].energy);

      // Then edit the event record
      const double scale = 1 + .03*sigma; // 3% energy scale syst.
      sr->sbn.truth.neutrino[0].energy *= scale;
    }
  };
  const ToyEnergyScaleSyst eSyst;

  // Make systematically shifted variants of the spectrum above
  Spectrum sEnergyUp(loaderBeam, axEnergy, kSelectionCut, SystShifts(&eSyst, +1));
  Spectrum sEnergyDn(loaderBeam, axEnergy, kSelectionCut, SystShifts(&eSyst, -1));

  class ToyNormSyst: public ISyst
  {
  public:
    ToyNormSyst() : ISyst("toyNorm", "Toy Norm Scale") {}

    void Shift(double sigma,
               Restorer& restore,
               caf::StandardRecord* sr,
               double& weight) const override
    {
      // A systematic can also reweight events, based on whatever criteria you
      // want.
      if(sr->sbn.truth.neutrino[0].energy > 1.5) weight *= 1+0.2*sigma;
      else weight *= 1+0.1*sigma;
    }
  };
  const ToyNormSyst nSyst;

  Spectrum sNormUp(loaderBeam, axEnergy, kSelectionCut, SystShifts(&nSyst, +1));
  Spectrum sNormDn(loaderBeam, axEnergy, kSelectionCut, SystShifts(&nSyst, -1));

  // Fill all the various shifted spectra
  loaderBeam.Go();

  // Plot what we have so far
  TCanvas* c1 = new TCanvas("c1");

  sEnergyDn.ToTH1(pot, kBlue)->Draw("hist");
  sEnergyUp.ToTH1(pot, kRed)->Draw("hist same");
  sEnergy.ToTH1(pot)->Draw("hist same");
  c1->SaveAs("demo4_plot1.pdf");

  // Now the normalisation systematic
  sNormUp.ToTH1(pot, kRed)->Draw("hist");
  sNormDn.ToTH1(pot, kBlue)->Draw("hist same");
  sEnergy.ToTH1(pot)->Draw("hist same");
  c1->SaveAs("demo4_plot2.pdf");
}
