#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Analysis/Calcs.h"
#include "CAFAna/Core/OscCalcSterileApprox.h"
#include "OscLib/func/OscCalculatorSterile.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "StandardRecord/Proxy/SRProxy.h"
#include "CAFAna/Core/Loaders.h"
#include "CAFAna/Prediction/PredictionInterp.h"
#include "CAFAna/Prediction/PredictionGenerator.h"
#include "TFile.h"
#include "CAFAna/Analysis/ExpInfo.h"

#include "toysysts.h"

// Random numbers
#include "TRandom3.h"

using namespace ana;

const std::string numuStr = "numu";
const std::string nueStr = "nue";

void make_state(const std::string anatype = numuStr)
{

  Loaders loaders, loaders2;
  if (anatype == numuStr) {
    //const std::string fDir = "/pnfs/sbnd/persistent/users/gputnam/numu_simulation_reweight/processed_2.a/";
    const std::string fDir = "/pnfs/sbnd/persistent/users/gputnam/numu_simulation_12_05_2018/processed_1.tempwgh/";
    const std::string fnameBeam = fDir + "output_SBNOsc_NumuSelection_Modern_SBND.root";
    const std::string fnameBeam2 = fDir + "output_SBNOsc_NumuSelection_Modern_Icarus.root";

    loaders.SetLoaderPath( fnameBeam, Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);
    loaders2.SetLoaderPath( fnameBeam2, Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);
  }
  else if (anatype == nueStr) {
    //std::cout << "Nue files not working right now!" << std::endl;
    //return;
    const std::string fDir = "/pnfs/sbn/persistent/users/dbarker/sbnoutput/";

    //BNB files contain nominal non-swap beam (so numubg, nuebg, NC)
    const std::string fnameBNB = fDir + "output_SBNOsc_NueSelection_Proposal_SBND_Numu.root";
    const std::string fnameBNB2 = fDir + "output_SBNOsc_NueSelection_Proposal_Icarus_NuMu.root";

    //Nue instrinsic only (to increase stats)
    const std::string fnameIntrinsic = fDir + "output_SBNOsc_NueSelection_Proposal_SBND_Int.root";
    const std::string fnameIntrinsic2 = fDir + "output_SBNOsc_NueSelection_Proposal_Icarus_Int.root";

    //Swap files are for signal
    const std::string fnameSwap = fDir + "output_SBNOsc_NueSelection_Proposal_SBND_Osc.root";
    const std::string fnameSwap2 = fDir + "output_SBNOsc_NueSelection_Proposal_Icarus_Osc.root";

    loaders.SetLoaderPath( fnameBNB,  Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);
    loaders.SetLoaderPath( fnameIntrinsic,  Loaders::kMC,   ana::kBeam, Loaders::kIntrinsic);
    loaders.SetLoaderPath( fnameSwap,  Loaders::kMC,   ana::kBeam, Loaders::kNueSwap);
    loaders2.SetLoaderPath( fnameBNB2, Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);
    loaders2.SetLoaderPath( fnameIntrinsic2,  Loaders::kMC,   ana::kBeam, Loaders::kIntrinsic);
    loaders2.SetLoaderPath( fnameSwap2,  Loaders::kMC,   ana::kBeam, Loaders::kNueSwap);

  }
  else {
    std::cout << "Unrecognized analysis - use numu or nue" << std::endl;
    return;
  }

  const double sbndPOT = kPOTnominal;
  const double icarusPOT = kPOTnominal;

  // Calculator (just use no oscillations)
  osc::NoOscillations* calc = new osc::NoOscillations;

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

  const Var kWeighthack([](const caf::SRProxy* sr)
                        {
			  if (sr->truth[0].neutrino.iscc && sr->truth[0].neutrino.pdg == 12) return 0.8;
			  if (sr->truth[0].neutrino.iscc && sr->truth[0].neutrino.pdg == 14) return 0.005;
			  if (sr->truth[0].neutrino.isnc) return 0.05;
			  return 1.0;
			});

  const Binning binsEnergy = Binning::Simple(30, 0, 3);
  const HistAxis axEnergy("Reconstructed energy (GeV)", binsEnergy, kRecoE);
  const HistAxis axTrueEnergy("True energy (GeV)", binsEnergy, kTrueE);

  // List all of the systematics we'll be using
  std::cout << "\nIncluding the following systematics:" << std::endl;
  for(const ISyst* s: allSysts) std::cout << s->ShortName() << "\t\t" << s->LatexName() << std::endl;
  std::cout << "\n" << std::endl;
  std::vector<const ISyst*> noSysts{};

  //Use true energy, no weights until we get new nue files
  NoExtrapGenerator gen(anatype == numuStr ? axEnergy : axTrueEnergy,
                        kNoCut,
			anatype == numuStr ? kWeight : kWeighthack);
  if (anatype == numuStr) {
    std::cout << "Using reco energy" << std::endl;
  }
  else {
    std::cout << "Using true energy" << std::endl;
  }

  PredictionInterp pred_nd(anatype == numuStr ? allSysts : noSysts, calc, gen, loaders);
  PredictionInterp pred_fd(anatype == numuStr ? allSysts : noSysts, calc, gen, loaders2);

  loaders.Go();
  loaders2.Go();

  std::cout << "Creating file " << ("cafe_state_smear_"+anatype+".root").c_str() << std::endl;

  TFile fout(("cafe_state_smear_"+anatype+".root").c_str(), "RECREATE");
  pred_nd.SaveTo(fout.mkdir(("pred_nd_"+anatype).c_str()));
  pred_fd.SaveTo(fout.mkdir(("pred_fd_"+anatype).c_str()));

}
