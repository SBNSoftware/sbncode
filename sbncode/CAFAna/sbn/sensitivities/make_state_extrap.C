#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "StandardRecord/Proxy/SRProxy.h"
#include "CAFAna/Core/Loaders.h"
#include "TFile.h"
#include "CAFAna/Analysis/ExpInfo.h"

#include "CAFAna/Prediction/PredictionSBNExtrap.h"

using namespace ana;

void make_state_extrap()
{

  Loaders loaders, loaders2;
  const std::string fDir = "/sbnd/data/users/bckhouse/processed_1.tempwgh/";
  const std::string fnameBeam = fDir + "output_SBNOsc_NumuSelection_Modern_SBND.root";
  const std::string fnameBeam2 = fDir + "output_SBNOsc_NumuSelection_Modern_Icarus.root";

  loaders.SetLoaderPath( fnameBeam, Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);

  // Also fake ND data
  loaders.SetLoaderPath( fnameBeam, Loaders::kData,   ana::kBeam, Loaders::kNonSwap);

  loaders2.SetLoaderPath( fnameBeam2, Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);

  const double sbndPOT = kPOTnominal;
  const double icarusPOT = kPOTnominal;

  const Var kRecoE([](const caf::SRProxy* sr)
                   {
                     return sr->reco[0].reco_energy;
                   });

  const Var kTrueE([](const caf::SRProxy* sr)
                        {
			  return sr->truth[0].neutrino.energy;
			});

  const Var kWeight([](const caf::SRProxy* sr)
		    {
		      return sr->reco[0].weight;
		    });

  const Var kWeightHack([](const caf::SRProxy* sr)
                        {
			  if (sr->truth[0].neutrino.iscc && sr->truth[0].neutrino.pdg == 12) return 0.8;
			  if (sr->truth[0].neutrino.iscc && sr->truth[0].neutrino.pdg == 14) return 0.005;
			  if (sr->truth[0].neutrino.isnc) return 0.05;
			  return 1.0;
			});

  const Binning binsEnergy = Binning::Simple(30, 0, 3);
  const HistAxis axEnergy("Reconstructed energy (GeV)", binsEnergy, kRecoE);
  const HistAxis axTrueEnergy("True energy (GeV)", binsEnergy, kTrueE);

  PredictionSBNExtrap pred(loaders, loaders2, axEnergy, kNoCut);

  PredictionNoExtrap predND(loaders, axEnergy, kNoCut);

  loaders.Go();
  loaders2.Go();

  std::cout << "Creating file " << "cafe_state_extrap.root" << std::endl;

  TFile fout("cafe_state_extrap.root", "RECREATE");
  pred.SaveTo(fout.mkdir("pred"));
  predND.SaveTo(fout.mkdir("predND"));
}
