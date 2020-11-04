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
#include "CAFAna/Systs/SBNFluxSysts.h"

#include "OscLib/IOscCalc.h"

#include "TFile.h"

#include <string>

using namespace ana;

const std::string numuStr = "numu";
const std::string nueStr = "nue";

void make_state_syst(const std::string anatype = numuStr)
{
  Loaders loaders_nd, loaders_fd, loaders_ub;

  if(anatype == numuStr) {
    const std::string dir = "/sbnd/data/users/jlarkin/workshop_samples/";
    const std::string fnameBeam_nd = dir + "output_SBNOsc_NumuSelection_Modern_SBND.flat.root";
    const std::string fnameBeam_fd = dir + "output_SBNOsc_NumuSelection_Modern_Icarus.flat.root";
    const std::string fnameBeam_ub = dir + "output_SBNOsc_NumuSelection_Modern_Uboone.flat.root";

    loaders_nd.SetLoaderPath(fnameBeam_nd, Loaders::kMC, ana::kBeam, Loaders::kNonSwap);
    loaders_fd.SetLoaderPath(fnameBeam_fd, Loaders::kMC, ana::kBeam, Loaders::kNonSwap);
    loaders_ub.SetLoaderPath(fnameBeam_ub, Loaders::kMC, ana::kBeam, Loaders::kNonSwap);
  }

  else if (anatype == nueStr) {
    const std::string dir = "/sbnd/data/users/bzamoran/Modern_6thNov2019/";

    //BNB files contain nominal non-swap beam (so numubg, nuebg, NC)
    const std::string fnameBNB_nd = dir + "output_SBNOsc_NueSelection_Modern_SBND_Numu.flat.root";
    const std::string fnameBNB_fd = dir + "output_SBNOsc_NueSelection_Modern_Icarus_Numu.flat.root";
    const std::string fnameBNB_ub = dir + "output_SBNOsc_NueSelection_Modern_Uboone_Numu.flat.root";

    //Nue instrinsic only (to increase stats)
    const std::string fnameInt_nd = dir + "output_SBNOsc_NueSelection_Modern_SBND_Int.flat.root";
    const std::string fnameInt_fd = dir + "output_SBNOsc_NueSelection_Modern_Icarus_Int.flat.root";
    const std::string fnameInt_ub = dir + "output_SBNOsc_NueSelection_Modern_Uboone_Int.flat.root";

    //Swap files are for signal
    const std::string fnameSwap_nd = dir + "output_SBNOsc_NueSelection_Modern_SBND_Osc.flat.root";
    const std::string fnameSwap_fd = dir + "output_SBNOsc_NueSelection_Modern_Icarus_Osc.flat.root";
    const std::string fnameSwap_ub = dir + "output_SBNOsc_NueSelection_Modern_Uboone_Osc.flat.root";

    // //Dirt (background)
    // const std::string fnameDirt_nd = dir + "output_SBNOsc_NueSelection_Modern_SBND_CosDirt.flat.root";
    // const std::string fnameDirt_fd = dir + "output_SBNOsc_NueSelection_Modern_Icarus_CosDirt.flat.root";
    // const std::string fnameDirt_ub = dir + "output_SBNOsc_NueSelection_Modern_Uboone_CosDirt.flat.root";

    loaders_nd.SetLoaderPath( fnameBNB_nd,   Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);
    loaders_nd.SetLoaderPath( fnameInt_nd,   Loaders::kMC,   ana::kBeam, Loaders::kIntrinsic);
    loaders_nd.SetLoaderPath( fnameSwap_nd,  Loaders::kMC,   ana::kBeam, Loaders::kNueSwap);

    loaders_fd.SetLoaderPath( fnameBNB_fd,   Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);
    loaders_fd.SetLoaderPath( fnameInt_fd,   Loaders::kMC,   ana::kBeam, Loaders::kIntrinsic);
    loaders_fd.SetLoaderPath( fnameSwap_fd,  Loaders::kMC,   ana::kBeam, Loaders::kNueSwap);

    loaders_ub.SetLoaderPath( fnameBNB_ub,   Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);
    loaders_ub.SetLoaderPath( fnameInt_ub,   Loaders::kMC,   ana::kBeam, Loaders::kIntrinsic);
    loaders_ub.SetLoaderPath( fnameSwap_ub,  Loaders::kMC,   ana::kBeam, Loaders::kNueSwap);
  }
  else {
    std::cout << "Unrecognised analysis - use numu or nue" << std::endl;
    return;
  }

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

  NoExtrapGenerator nom_gen(axEnergy, kOneTrue, kWeight);

  std::vector<const ISyst*> systs = GetSBNWeightSysts();
  for(const ISyst* s: GetSBNFluxHadronSysts(30)) systs.push_back(s);

  osc::NoOscillations calc;

  PredictionInterp pred_nd(systs, &calc, nom_gen, loaders_nd);
  PredictionInterp pred_fd(systs, &calc, nom_gen, loaders_fd);
  PredictionInterp pred_ub(systs, &calc, nom_gen, loaders_ub);

  loaders_nd.Go();
  loaders_fd.Go();
  loaders_ub.Go();

  std::cout << "Creating file " << ("cafe_state_syst_"+anatype+".root").c_str() << std::endl;

  TFile fout(("cafe_state_syst_"+anatype+".root").c_str(), "RECREATE");

  pred_nd.SaveTo(fout.mkdir("pred_nd"));
  pred_fd.SaveTo(fout.mkdir("pred_fd"));
  pred_ub.SaveTo(fout.mkdir("pred_ub"));
}



