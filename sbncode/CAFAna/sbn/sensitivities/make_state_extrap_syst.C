#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Prediction/PredictionInterp.h"
#include "StandardRecord/Proxy/SRProxy.h"
#include "CAFAna/Core/Loaders.h"
#include "TFile.h"
#include "CAFAna/Analysis/ExpInfo.h"

#include "CAFAna/Systs/SBNWeightSysts.h"
#include "CAFAna/Systs/SBNFluxSysts.h"

#include "OscLib/IOscCalc.h"

#include "CAFAna/Prediction/PredictionSBNExtrap.h"

using namespace ana;

void make_state_extrap_syst()
{
  Loaders loaders_nd, loaders_fd;
  const std::string fDir = "/sbnd/data/users/jlarkin/workshop_samples/";
  const std::string fnameBeam_nd = "my.sbnd.flat.lz4.1.root"; // fDir + "output_SBNOsc_NumuSelection_Modern_SBND.flat.root";
  const std::string fnameBeam_fd = "my.icarus.flat.lz4.1.root"; // fDir + "output_SBNOsc_NumuSelection_Modern_Icarus.flat.root";

  loaders_nd.SetLoaderPath( fnameBeam_nd, Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);

  // Also fake ND data
  loaders_nd.SetLoaderPath( fnameBeam_nd, Loaders::kData,   ana::kBeam, Loaders::kNonSwap);

  loaders_fd.SetLoaderPath( fnameBeam_fd, Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);

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


  std::vector<const ISyst*> systs = GetSBNWeightSysts();
  for(const ISyst* s: GetSBNFluxHadronSysts(30)) systs.push_back(s);

  osc::NoOscillations calc;

  NoExtrapGenerator nogen(axEnergy, kOneTrue, kWeight);
  SBNExtrapGenerator extrapgen(loaders_nd, axEnergy, kOneTrue, kWeight,
                               /*for data*/ kNoShift, kWeight);


  PredictionSBNExtrap pred(loaders_nd, loaders_fd, axEnergy, kOneTrue, kNoShift, kWeight);

  PredictionInterp predInterp(systs, &calc, extrapgen, loaders_fd);

  PredictionNoExtrap predND(loaders_nd, axEnergy, kOneTrue, kNoShift, kWeight);

  PredictionInterp predNDInterp(systs, &calc, nogen, loaders_nd);

  loaders_fd.Go();
  loaders_nd.Go();

  std::cout << "Creating file " << "cafe_state_extrap.root" << std::endl;

  TFile fout("cafe_state_extrap_syst.root", "RECREATE");
  pred.SaveTo(fout.mkdir("pred"));
  predInterp.SaveTo(fout.mkdir("pred_syst"));
  predND.SaveTo(fout.mkdir("predND"));
  predNDInterp.SaveTo(fout.mkdir("predND_syst"));
}
