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

void make_state_smear()
{
  const std::string fDir = "/pnfs/sbnd/persistent/users/gputnam/numu_simulation_reweight/processed_2.a/";
  // First SBND data
  const std::string fnameBeam = fDir + "output_SBNOsc_NumuSelection_Modern_SBND.root";
  // And now add Icarus data
  const std::string fnameBeam2 = fDir + "output_SBNOsc_NumuSelection_Modern_Icarus.root";

  // Source of events
  Loaders loaders;
  loaders.SetLoaderPath( fnameBeam,  caf::kFARDET,  Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);
  //loaders.SetLoaderPath( fnameSwap,  caf::kFARDET,  Loaders::kMC,   ana::kBeam, Loaders::kFluxSwap);
  Loaders loaders2;
  loaders2.SetLoaderPath( fnameBeam2,  caf::kFARDET,  Loaders::kMC,   ana::kBeam, Loaders::kNonSwap);

  const double sbndPOT = POTnominal;
  const double icarusPOT = POTnominal;

  // Calculator
  osc::OscCalculatorSterile* calc = DefaultSterileCalc(4);
  calc->SetL(BaselineSBND);
  osc::OscCalculatorSterile* calc2 = DefaultSterileCalc(4);
  calc2->SetL(BaselineIcarus);

  // This is probably too simplistic. Maybe res/sqrt(E)?
  const Var kSmearedE([](const caf::SRProxy* sr)
                        {
			  return sr->reco[0].reco_energy;
			});

  const Var kWeight([](const caf::SRProxy* sr)
                        {
			  return sr->reco[0].weight;
			});

  const Binning binsEnergy = Binning::Simple(30, 0, 3);
  const HistAxis axEnergy("Fake reconstructed energy (GeV)", binsEnergy, kSmearedE);

  // List all of the systematics we'll be using
  std::cout << "\nIncluding the following systematics:" << std::endl;
  for(const ISyst* s: allSysts) std::cout << s->ShortName() << "\t\t" << s->LatexName() << std::endl;
  std::cout << "\n" << std::endl;

  NoExtrapGenerator gen(axEnergy, kNoCut, kWeight);

  PredictionInterp pred_nd_numu(allSysts, calc, gen, loaders);
  PredictionInterp pred_fd_numu(allSysts, calc2, gen, loaders2);

  loaders.Go();
  loaders2.Go();

  const char* stateFname = "cafe_state_smear.root";

  TFile fout(stateFname, "RECREATE");
  pred_nd_numu.SaveTo(fout.mkdir("pred_nd_numu"));
  pred_fd_numu.SaveTo(fout.mkdir("pred_fd_numu"));

}
