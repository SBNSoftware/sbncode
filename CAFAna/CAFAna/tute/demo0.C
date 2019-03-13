// Make a simple spectrum plot
// cafe demo0.C

#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"

// #include "CAFAna/Cuts/TruthCuts.h"

#include "StandardRecord/StandardRecord.h"

#include "TCanvas.h"
#include "TH1.h"
#include "TPad.h"

using namespace ana;

void demo0()
{
  // Environment variables and wildcards work. As do SAM datasets.
  const std::string fname = "/sbnd/app/users/bzamoran/sbncode-v07_11_00/output_largesample_nu_ExampleAnalysis_ExampleSelection.root";

  // Source of events
  SpectrumLoader loader(fname);

  // A Var is a little snippet of code that takes a record representing the
  // event record and returns a single number to plot.
  const Var kTruthEnergy({},
                        [](const caf::StandardRecord* sr)
                        {
                          return sr->sbn.truth.neutrino[0].energy;
                        });

  // For such a simple variable you can use a shortcut like this
  const Var kTruthY = SIMPLEVAR(sbn.truth.neutrino[0].inelasticityY);

  // Define a spectrum, ie a histogram with associated POT information
  const Binning binsEnergy = Binning::Simple(50, 0, 5);
  const HistAxis axEnergy("True energy (GeV)", binsEnergy, kTruthEnergy);
  // kIsNumuCC here is a "Cut". Same as a Var but returning a boolean. In this
  // case, we're only keeping events that are truly numu CC interactions.
  const Var kIsCC = SIMPLEVAR(sbn.truth.neutrino[0].iscc);
  Spectrum sEnergy(loader, axEnergy, kIsCC != 0);

  Spectrum sEnergyNC(loader, axEnergy, kIsCC == 0);

  // And another
  const Binning binsY = Binning::Simple(50, 0, 1);
  const HistAxis axY("Y", binsY, kTruthY);
  Spectrum sY(loader, axY, kIsCC != 0);
  Spectrum sYNC(loader, axY, kIsCC == 0);

  // This is the call that actually fills in those spectra
  loader.Go();

  // Fake POT: we need to sort this out in the files first
  const double pot = 6.e20;

  // For plotting purposes we can convert to TH1s
  TCanvas* c1 = new TCanvas("c1");
  sEnergy.ToTH1(pot)->Draw("hist");
  sEnergyNC.ToTH1(pot, kBlue)->Draw("hist same");
  c1->SaveAs("demo0_plot1.pdf");

  sY.ToTH1(pot)->Draw("hist");
  sYNC.ToTH1(pot, kBlue)->Draw("hist same");
  gPad->SetLogy();
  c1->SaveAs("demo0_plot2.pdf");
}
