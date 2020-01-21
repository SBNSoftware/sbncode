// Making multiple plots at once
// cafe demo0b.C

// Includes from demo0.C
#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Analysis/ExpInfo.h"
#include "StandardRecord/Proxy/SRProxy.h"

// New includes
#include "CAFAna/Core/Ratio.h"
#include "CAFAna/Cuts/TruthCuts.h" // for kIsNC

using namespace ana;

#include "TCanvas.h"
#include "TH1.h"
#include "TPad.h"

void demo0b()
{
  // We will make spectra from two different files (for two different
  // detectors)
  const std::string fnameSBND = "/sbnd/data/users/bckhouse/sample_2.1_fitters/output_SBNOsc_NumuSelection_Proposal_SBND.flat.root";
  const std::string fnameIcarus = "/sbnd/data/users/bckhouse/sample_2.1_fitters/output_SBNOsc_NumuSelection_Proposal_Icarus.flat.root";

  SpectrumLoader loaderSBND(fnameSBND);
  SpectrumLoader loaderIcarus(fnameIcarus);

  // We will also use two different Vars. For such simple ones we can use this
  // shortcut macro.
  const Var kTruthEnergy = SIMPLEVAR(truth[0].neutrino.energy);
  const Var kTruthY = SIMPLEVAR(truth[0].neutrino.inelasticityY);

  // Define axes for the spectra we'll make
  const HistAxis axEnergy("True energy (GeV)", Binning::Simple(50, 0, 5), kTruthEnergy);
  const HistAxis axY("Y", Binning::Simple(50, 0, 1), kTruthY);

  // A Cut is the same as a Var but returns a boolean. Here we're selecting
  // truly numu CC interactions.
  const Cut kIsNumuCC([](const caf::SRProxy* sr)
                      {
                        return (sr->truth[0].neutrino.iscc &&
                                abs(sr->truth[0].neutrino.pdg) == 14);
                      });

  // Can also make combinations of existing Cuts and Vars in the obvious
  // ways. kIsNC here comes from Cuts/TruthCuts.h. Really we should have just
  // used kIsNumuCC from there.
  const Var kPdg = SIMPLEVAR(truth[0].neutrino.pdg);
  const Cut kIsNumuCC2 = (kPdg == 14 || kPdg == -14) && !kIsNC;

  Spectrum sEnergySig(loaderSBND, axEnergy, kIsNumuCC);
  Spectrum sEnergyBkg(loaderSBND, axEnergy, !kIsNumuCC);

  Spectrum sEnergySigIc(loaderIcarus, axEnergy, kIsNumuCC);
  Spectrum sEnergyBkgIc(loaderIcarus, axEnergy, !kIsNumuCC);

  Spectrum sYSig(loaderSBND, axY, kIsNumuCC);
  Spectrum sYBkg(loaderSBND, axY, !kIsNumuCC);

  // Fill in the 4 ND spectra
  loaderSBND.Go();
  // And the 2 FD
  loaderIcarus.Go();

  sEnergySig.ToTH1(kPOTnominal)->Draw("hist");
  sEnergyBkg.ToTH1(kPOTnominal, kBlue)->Draw("hist same");

  new TCanvas;

  sYSig.ToTH1(kPOTnominal)->Draw("hist");
  sYBkg.ToTH1(kPOTnominal, kBlue)->Draw("hist same");

  new TCanvas;

  TH1* h = sEnergySigIc.ToTH1(kPOTnominal);
  h->SetTitle("Icarus");
  h->Draw("hist");
  sEnergyBkgIc.ToTH1(kPOTnominal, kBlue)->Draw("hist same");

  new TCanvas;

  // And just for fun...
  // Note no need to mention POT anywhere when forming a ratio
  Ratio r(sEnergySigIc, sEnergySig);
  r.ToTH1()->Draw();
}
