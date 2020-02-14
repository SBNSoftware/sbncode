// Oscillations using Prediction class
// cafe demo1b.C

#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/OscCalcSterileApprox.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Analysis/ExpInfo.h"

using namespace ana;

#include "StandardRecord/Proxy/SRProxy.h"

#include "TH1.h"

// New includes for this macro
#include "CAFAna/Prediction/PredictionNoExtrap.h"

void demo1b()
{
  // Repeated from previous macros
  const std::string fnameIcarus = "/sbnd/data/users/bckhouse/sample_2.1_fitters/output_SBNOsc_NumuSelection_Proposal_Icarus.flat.root";
  SpectrumLoader loaderIcarus(fnameIcarus);
  const Var kRecoE = SIMPLEVAR(reco.reco_energy);
  const Var kWeight = SIMPLEVAR(reco.weight);
  const HistAxis axEnergy("Reconstructed energy (GeV)", Binning::Simple(50, 0, 5), kRecoE);
  OscCalcSterileApprox calc;
  calc.SetDmsq(1);
  calc.SetSinSq2ThetaMuMu(0.2);
  calc.SetSinSq2ThetaMuE(0);
  calc.SetL(kBaselineIcarus);

  // A Prediction is an objects holding a variety of "OscillatableSpectrum"
  // objects, one for each original and final flavour combination.
  PredictionNoExtrap pred(loaderIcarus, kNullLoader, kNullLoader, kNullLoader,
                          axEnergy, kNoCut, kNoShift, kWeight);

  // This call will fill all of the constituent parts of the prediction
  loaderIcarus.Go();

  // We can extract a total prediction unoscillated
  const Spectrum sUnosc = pred.PredictUnoscillated();
  // Or oscillated
  const Spectrum sOsc = pred.Predict(&calc);

  // And we can break things down by flavour
  const Spectrum sOscNC = pred.PredictComponent(&calc,
                                                Flavors::kAll,
                                                Current::kNC,
                                                Sign::kBoth);

  sUnosc.ToTH1(kPOTnominal)->Draw("hist");
  sOsc.ToTH1(kPOTnominal, kRed)->Draw("hist same");
  sOscNC.ToTH1(kPOTnominal, kBlue)->Draw("hist same");
}
