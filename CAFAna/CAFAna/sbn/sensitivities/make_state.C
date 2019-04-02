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
    //kFARDET is NOvA residual
    loaders.SetLoaderPath( fnameBeam, Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);
    loaders2.SetLoaderPath( fnameBeam2, Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);
  }
  else if (anatype == nueStr) {
    std::cout << "Nue files not working right now!" << std::endl;
    return;
    const std::string fDir = "/pnfs/sbn/persistent/users/dbarker/sbnoutput/";

    //BNB files contain nominal non-swap beam (so numubg, nuebg, NC)
    const std::string fnameBNB = fDir + "output_SBNOsc_NueSelection_Proposal_SBND.root";
    const std::string fnameBNB2 = fDir + "output_SBNOsc_NueSelection_Proposal_Icarus.root";

    //Non-swap files with nue instrinsic only
    //const std::string fnameIntrinsic = fDir + "?";
    //const std::string fnameIntrinsic = fDir + "?";

    //Swap files are for signal
    //const std::string fnameSwap = fDir + "?";
    //const std::string fnameSwap2 = fDir + "?";

    loaders.SetLoaderPath( fnameBNB,  Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);
    //loaders.SetLoaderPath( fnameIntrinsic,  Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);
    //loaders.SetLoaderPath( fnameSwap,  Loaders::kMC,   ana::kBeam, Loaders::kNueSwap);
    loaders2.SetLoaderPath( fnameBNB2, Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);
    //loaders2.SetLoaderPath( fnameIntrinsic2,  Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);
    //loaders2.SetLoaderPath( fnameSwap2,  Loaders::kMC,   ana::kBeam, Loaders::kNueSwap);

  }
  else {
    std::cout << "Unrecognized analysis - use numu or nue" << std::endl;
    return;
  }

  const double sbndPOT = kPOTnominal;
  const double icarusPOT = kPOTnominal;

  // Calculator
  OscCalcSterileApproxAdjustable* calc = DefaultSterileApproxCalc();
  calc->SetL(kBaselineSBND);
  OscCalcSterileApproxAdjustable* calc2 = DefaultSterileApproxCalc();
  calc2->SetL(kBaselineIcarus);

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

  const Binning binsEnergy = Binning::Simple(30, 0, 3);
  const HistAxis axEnergy("Reconstructed energy (GeV)", binsEnergy, kRecoE);
  const HistAxis axTrueEnergy("True energy (GeV)", binsEnergy, kTrueE);

  // List all of the systematics we'll be using
  std::cout << "\nIncluding the following systematics:" << std::endl;
  for(const ISyst* s: allSysts) std::cout << s->ShortName() << "\t\t" << s->LatexName() << std::endl;
  std::cout << "\n" << std::endl;

  //Use true energy, no weights until we get new nue files
  NoExtrapGenerator gen(axTrueEnergy, kNoCut, kUnweighted);
  if (anatype == numuStr) {
    NoExtrapGenerator gen(axEnergy, kNoCut, kWeight);
  }

  PredictionInterp pred_nd(allSysts, calc, gen, loaders);
  PredictionInterp pred_fd(allSysts, calc2, gen, loaders2);

  loaders.Go();
  loaders2.Go();

  std::cout << "Creating file " << ("cafe_state_smear_"+anatype+".root").c_str() << std::endl;

  TFile fout(("cafe_state_smear_"+anatype+".root").c_str(), "RECREATE");
  pred_nd.SaveTo(fout.mkdir(("pred_nd_"+anatype).c_str()));
  pred_fd.SaveTo(fout.mkdir(("pred_fd_"+anatype).c_str()));

}
