// Oscillations
// cafe demo1.C

#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Analysis/ExpInfo.h"

// New includes for this macro
#include "CAFAna/Core/OscillatableSpectrum.h"
#include "CAFAna/Core/OscCalcSterileApprox.h"

using namespace ana;

#include "StandardRecord/Proxy/SRProxy.h"

#include "TCanvas.h"
#include "TH2.h"

void demo1()
{
  // Repeated from previous macros
  const std::string fnameSBND = "/sbnd/data/users/bckhouse/sample_2.1_fitters/output_SBNOsc_NumuSelection_Proposal_SBND.flat.root";
  const std::string fnameIcarus = "/sbnd/data/users/bckhouse/sample_2.1_fitters/output_SBNOsc_NumuSelection_Proposal_Icarus.flat.root";
  SpectrumLoader loaderSBND(fnameSBND);
  SpectrumLoader loaderIcarus(fnameIcarus);

  // Mocked up reconstruction
  const Var kRecoE = SIMPLEVAR(reco.reco_energy);
  // Mock selection efficiencies etc
  const Var kWeight = SIMPLEVAR(reco.weight);

  const HistAxis axEnergy("Reconstructed energy (GeV)", Binning::Simple(50, 0, 5), kRecoE);

  // These are 2D reco-vs-true distributions
  OscillatableSpectrum sOscND(loaderSBND,   axEnergy, kIsNumuCC, kNoShift, kWeight);
  OscillatableSpectrum sOscFD(loaderIcarus, axEnergy, kIsNumuCC, kNoShift, kWeight);

  loaderSBND.Go();
  loaderIcarus.Go();

  sOscFD.ToTH2(kPOTnominal)->Draw("colz");
  gPad->SetLogy();

  // Full three/four-flavour calculators are available, but for most cases we
  // want to use this simple two-flavour model, which also takes account of
  // rapid oscillations.
  OscCalcSterileApprox calc;

  calc.SetDmsq(1); // oscillations at 1eV^2

  calc.SetSinSq2ThetaMuMu(0.2); // numu disappearance
  calc.SetSinSq2ThetaMuE(0);

  // NB: the baseline is a property of the oscillation calculator! You will
  // need to set it appropriately before calls involving each detector.
  calc.SetL(kBaselineSBND);

  new TCanvas;

  // This is just the projection onto the reco axis
  const Spectrum sUnosc = sOscND.Unoscillated();
  // And this is the same projection, but applying oscillation weights (here numu survival)
  const Spectrum sOsc = sOscND.Oscillated(&calc, 14, 14);

  sUnosc.ToTH1(kPOTnominal)->Draw("hist");
  sOsc.ToTH1(kPOTnominal, kRed)->Draw("hist same");

  // Have to set this before making the corresponding Icarus plot
  calc.SetL(kBaselineIcarus);

  new TCanvas;

  sOscFD.Unoscillated().ToTH1(kPOTnominal)->Draw("hist");
  sOscFD.Oscillated(&calc, 14, 14).ToTH1(kPOTnominal, kRed)->Draw("hist same");
}
