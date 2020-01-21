// Make a simple spectrum plot
// cafe demo0.C

#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"

#include "CAFAna/Analysis/ExpInfo.h" // for kPOTnominal

using namespace ana;

#include "StandardRecord/Proxy/SRProxy.h"

#include "TH1.h"

void demo0()
{
  // Converted from the original sbnanalysis files here
  // /sbnd/data/users/gputnam/NuMu/outputs/sample_2.1_fitters/
  //
  // Environment variables and wildcards work as arguments to
  // SpectrumLoader. As do SAM datasets.
  const std::string fname = "/sbnd/data/users/bckhouse/sample_2.1_fitters/output_SBNOsc_NumuSelection_Proposal_SBND.flat.root";

  // Source of events
  SpectrumLoader loader(fname);

  // A Var is a little snippet of code that takes a "proxy" representing the
  // event record and returns a single number to plot
  const Var kTruthEnergy([](const caf::SRProxy* sr)
                         {
                           return sr->truth[0].neutrino.energy;
                         });

  // Define a spectrum, ie a histogram with associated POT information
  const Binning binsEnergy = Binning::Simple(50, 0, 5);
  const HistAxis axEnergy("True energy (GeV)", binsEnergy, kTruthEnergy);
  Spectrum sEnergy(loader, axEnergy, kNoCut);

  // This is the call that actually fills in the spectrum
  loader.Go();

  // From Analysis/ExpInfo.h
  const double pot = kPOTnominal; // currently 6.6e20

  // For plotting purposes we can convert to a TH1
  sEnergy.ToTH1(pot)->Draw("hist");
}
