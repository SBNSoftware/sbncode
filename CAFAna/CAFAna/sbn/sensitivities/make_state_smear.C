#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Analysis/Calcs.h"
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

void make_state_smear(const std::string anatype = numuStr)
{

  Loaders loaders, loaders2;
  if (anatype == numuStr) {
    const std::string fDir = "/pnfs/sbnd/persistent/users/gputnam/numu_simulation_reweight/processed_2.a/";
    const std::string fnameBeam = fDir + "output_SBNOsc_NumuSelection_Modern_SBND.root";
    const std::string fnameBeam2 = fDir + "output_SBNOsc_NumuSelection_Modern_Icarus.root";
    //kFARDET is NOvA residual
    loaders.SetLoaderPath( fnameBeam,  caf::kFARDET,  Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);
    loaders2.SetLoaderPath( fnameBeam2,  caf::kFARDET,  Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);
  }
  else if (anatype == nueStr) {
    const std::string fDir = "/pnfs/sbn/persistent/users/dbarker/sbnoutput/";
    const std::string fnameSwap = fDir + "output_SBNOsc_NueSelection_Proposal_SBND.root";
    const std::string fnameSwap2 = fDir + "output_SBNOsc_NueSelection_Proposal_Icarus.root";
    loaders.SetLoaderPath( fnameSwap,  caf::kFARDET,  Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);
    loaders2.SetLoaderPath( fnameSwap2,  caf::kFARDET,  Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);
  }
  else {
    std::cout << "Unrecognized analysis - use numu or nue" << std::endl;
    return;
  }

  const double sbndPOT = kPOTnominal;
  const double icarusPOT = kPOTnominal;

  // Calculator
  osc::OscCalculatorSterile* calc = DefaultSterileCalc(4);
  calc->SetL(kBaselineSBND);
  osc::OscCalculatorSterile* calc2 = DefaultSterileCalc(4);
  calc2->SetL(kBaselineIcarus);

  // This is probably too simplistic. Maybe res/sqrt(E)?
  const Var kSmearedE([](const caf::SRProxy* sr)
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
  const HistAxis axEnergy("Reconstructed energy (GeV)", binsEnergy, kSmearedE);
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
