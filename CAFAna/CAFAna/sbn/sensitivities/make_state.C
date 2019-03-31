#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "StandardRecord/Proxy/SRProxy.h"
#include "TFile.h"

using namespace ana;

void make_state()
{
  const std::string fDir = "/pnfs/sbnd/persistent/users/gputnam/numu_simulation_reweight/processed_2.a/";

  const std::string fnameBeam = fDir + "output_SBNOsc_NumuSelection_Modern_SBND.root";

  // Source of events
  SpectrumLoader loaderBeam(fnameBeam);
  // SpectrumLoader loaderSwap(fnameSwap);

  // And now add Icarus data
  const std::string fnameBeam2 = fDir + "output_SBNOsc_NumuSelection_Modern_Icarus.root";

  // Source of events
  SpectrumLoader loaderBeam2(fnameBeam2);
  // SpectrumLoader loaderSwap2(fnameSwap2);

  const double sbndPOT = 6.6e20;
  const double icarusPOT = 6.6e20;

  //const Var kTrueE = SIMPLEVAR(truth.neutrino[0].energy);

  const Var kTrueE([](const caf::SRProxy* sr)
                         {
                           return sr->truth[0].neutrino.energy;
                         });

  const Binning binsEnergy = Binning::Simple(30, 0, 3);
  const HistAxis axEnergy("True energy (GeV)", binsEnergy, kTrueE);

  PredictionNoExtrap pred_nd_numu(loaderBeam, kNullLoader, kNullLoader,
                          axEnergy, kNoCut);

  PredictionNoExtrap pred_fd_numu(loaderBeam2, kNullLoader, kNullLoader,
                          axEnergy, kNoCut);

  loaderBeam.Go();
  // loaderSwap.Go();

  loaderBeam2.Go();
  // loaderSwap2.Go();

  const char* stateFname = "cafe_state.root";

  TFile fout(stateFname, "RECREATE");
  pred_nd_numu.SaveTo(fout.mkdir("pred_nd_numu"));
  pred_fd_numu.SaveTo(fout.mkdir("pred_fd_numu"));

}
