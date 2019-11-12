#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Analysis/Calcs.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "StandardRecord/Proxy/SRProxy.h"
#include "CAFAna/Core/Loaders.h"
#include "CAFAna/Prediction/PredictionInterp.h"
#include "CAFAna/Prediction/PredictionGenerator.h"

#include "CAFAna/Systs/SBNWeightSysts.h"

#include "OscLib/IOscCalculator.h"

#include "TFile.h"

#include <string>

using namespace ana;

void make_state_syst()
{
  const std::string dir = "/sbnd/data/users/jlarkin/workshop_samples/";
  const std::string fnameBeam_nd = dir + "output_SBNOsc_NumuSelection_Modern_SBND.flat.root";
  const std::string fnameBeam_fd = dir + "output_SBNOsc_NumuSelection_Modern_Icarus.flat.root";
  const std::string fnameBeam_ub = dir + "output_SBNOsc_NumuSelection_Modern_Uboone.flat.root";

  Loaders loaders_nd, loaders_fd, loaders_ub;

  loaders_nd.SetLoaderPath(fnameBeam_nd, Loaders::kMC, ana::kBeam, Loaders::kNonSwap);
  loaders_fd.SetLoaderPath(fnameBeam_fd, Loaders::kMC, ana::kBeam, Loaders::kNonSwap);
  loaders_ub.SetLoaderPath(fnameBeam_ub, Loaders::kMC, ana::kBeam, Loaders::kNonSwap);

  const Var kRecoE = SIMPLEVAR(reco.reco_energy);
  const Var kWeight = SIMPLEVAR(reco.weight);

  const vector<double> binEdges = {0.2, 0.3, 0.4, 0.45, 0.5,
                           0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0,
                           1.25, 1.5, 2.0, 2.5, 3.0};
  const Binning binsEnergy = Binning::Custom(binEdges);
  const HistAxis axEnergy("Reconstructed energy (GeV)", binsEnergy, kRecoE);

  NoExtrapGenerator nom_gen(axEnergy, kNoCut, kWeight);

  const std::vector<const ISyst*>& systs = GetSBNWeightSysts();

  osc::NoOscillations calc;

  PredictionInterp pred_nd(systs, &calc, nom_gen, loaders_nd);
  PredictionInterp pred_fd(systs, &calc, nom_gen, loaders_fd);
  PredictionInterp pred_ub(systs, &calc, nom_gen, loaders_ub);

  loaders_nd.Go();
  loaders_fd.Go();
  loaders_ub.Go();

  std::cout << "Creating file " << "cafe_state_syst_numu.root" << std::endl;

  TFile fout("cafe_state_syst_numu.root", "RECREATE");

  pred_nd.SaveTo(fout.mkdir("pred_nd"));
  pred_fd.SaveTo(fout.mkdir("pred_fd"));
  pred_ub.SaveTo(fout.mkdir("pred_ub"));
}



