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
  const std::string fDir = "/sbnd/data/users/jlarkin/workshop_samples/";
  const std::string fnameBeam = fDir + "output_SBNOsc_NumuSelection_Modern_SBND.flat.root";
  const std::string fnameBeam2 = fDir + "output_SBNOsc_NumuSelection_Modern_Icarus.flat.root";

  loaders.SetLoaderPath( fnameBeam, Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);

  // Also fake ND data
  loaders.SetLoaderPath( fnameBeam, Loaders::kData,   ana::kBeam, Loaders::kNonSwap);

  loaders2.SetLoaderPath( fnameBeam2, Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);

  const double sbndPOT = kPOTnominal;
  const double icarusPOT = kPOTnominal;

  const Var kRecoE = SIMPLEVAR(reco.reco_energy);
  const Var kWeight = SIMPLEVAR(reco.weight);
  const Cut kOneTrue([](const caf::SRProxy* sr)
         {
           return (sr->truth.size() == 1);
         });

  const vector<double> binEdges = {0.2, 0.3, 0.4, 0.45, 0.5,
                           0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0,
                           1.25, 1.5, 2.0, 2.5, 3.0};
  const Binning binsEnergy = Binning::Custom(binEdges);
  const HistAxis axEnergy("Reconstructed energy (GeV)", binsEnergy, kRecoE);

  PredictionSBNExtrap pred(loaders, loaders2, axEnergy, kOneTrue, kNoShift, kWeight);

  PredictionNoExtrap predND(loaders, axEnergy, kOneTrue, kNoShift, kWeight);

  loaders.Go();
  loaders2.Go();

  std::cout << "Creating file " << "cafe_state_extrap.root" << std::endl;

  TFile fout("cafe_state_extrap.root", "RECREATE");
  pred.SaveTo(fout.mkdir("pred"));
  predND.SaveTo(fout.mkdir("predND"));
}
