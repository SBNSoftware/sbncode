// Exercise the fitter
// cafe demo2.C

#include "CAFAna/Core/SpectrumLoader.h"
#include "CAFAna/Core/Spectrum.h"
#include "CAFAna/Core/Binning.h"
#include "CAFAna/Core/Var.h"
#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Core/OscCalcSterileApprox.h"
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Analysis/ExpInfo.h"

// New includes required
#include "CAFAna/Experiment/SingleSampleExperiment.h"
#include "CAFAna/Analysis/Fit.h"
#include "CAFAna/Analysis/FitAxis.h"
#include "CAFAna/Analysis/Surface.h"
#include "CAFAna/Vars/FitVarsSterileApprox.h"
#include "CAFAna/Experiment/MultiExperimentSBN.h"

using namespace ana;

#include "StandardRecord/Proxy/SRProxy.h"

#include "TCanvas.h"
#include "TH1.h"

void demo2(int step = 99)
{
  // Repeated from previous macros
  const std::string fnameIcarus = "/sbnd/data/users/bckhouse/sample_2.1_fitters/output_SBNOsc_NumuSelection_Proposal_Icarus.flat.root";
  const std::string fnameSBND = "/sbnd/data/users/bckhouse/sample_2.1_fitters/output_SBNOsc_NumuSelection_Proposal_SBND.flat.root";
  SpectrumLoader loaderIcarus(fnameIcarus);
  SpectrumLoader loaderSBND(fnameSBND);
  const Var kRecoE = SIMPLEVAR(reco.reco_energy);
  const Var kWeight = SIMPLEVAR(reco.weight);
  const HistAxis axEnergy("Reconstructed energy (GeV)", Binning::Simple(50, 0, 5), kRecoE);
  PredictionNoExtrap predSBND(loaderSBND, kNullLoader, kNullLoader, kNullLoader,
                              axEnergy, kNoCut, kNoShift, kWeight);
  PredictionNoExtrap predIcarus(loaderIcarus, kNullLoader, kNullLoader, kNullLoader,
                                axEnergy, kNoCut, kNoShift, kWeight);

  loaderSBND.Go();
  loaderIcarus.Go();


  // To make a fit we need to have a "data" spectrum to compare to our MC
  // Prediction object. "FakeData" here is an Asimov dataset. "MockData" would
  // include random poisson fluctuations.
  const Spectrum dataSBND = predSBND.PredictUnoscillated().FakeData(kPOTnominal);
  const Spectrum dataIcarus = predIcarus.PredictUnoscillated().FakeData(kPOTnominal);

  // An Experiment object is something that can turn oscillation parameters
  // into a chisq, in this case by comparing a Prediction and a data Spectrum
  SingleSampleExperiment exptSBND(&predSBND, dataSBND);
  SingleSampleExperiment exptIcarus(&predIcarus, dataIcarus);

  // A bit of faff setting this up this time because we need a version coerced
  // into satisfying the standard fitter interface.
  OscCalcSterileApproxAdjustable calc;
  calc.calc.SetDmsq(1);
  calc.calc.SetSinSq2ThetaMuMu(0.2); // these two don't really matter, the fitter will vary them.
  calc.calc.SetSinSq2ThetaMuE(0); // We'll be holding this one fixed at zero though

  // Let's make an SBND contour first
  calc.SetL(kBaselineSBND);

  // A fitter finds the minimum chisq using MINUIT by varying the list of
  // parameters given. These are FitVars from Vars/FitVars.h. They can contain
  // snippets of code to convert from the underlying angles etc to whatever
  // function you want to fit.
  Fitter fit(&exptSBND, {&kFitDmSqSterile, &kFitSinSq2ThetaMuMu});
  const double best_chisq = fit.Fit(&calc);

  if(step < 2) return;

  // 'calc' has been updated in place, so we could extract the best fit
  // oscillation parameters from it here

  //Define fit axes. 'true' here indicates a log scale
  const FitAxis kAxSinSq2ThetaMuMu(&kFitSinSq2ThetaMuMu, 30, 1e-3, 1, true);
  const FitAxis kAxDmSq(&kFitDmSqSterile, 30, 2e-2, 1e2, true);

  // A Surface is a map of chisq values in a 2D space
  Surface surfSBND(&exptSBND, &calc, kAxSinSq2ThetaMuMu, kAxDmSq);

  // Inspect the chisq map directly
  surfSBND.Draw();

  if(step < 3) return;

  // To draw a confidence interval we need a suitable critical value
  // surface. This generality is so that we can support Feldman-Cousins
  // corrections. For now we just want a constant critical value. Please check
  // carefully what significance, how many DoF, and one/two-sidedness you want.
  TH2* critSurf = Gaussian3Sigma1D1Sided(surfSBND);

  surfSBND.DrawContour(critSurf, kSolid, kRed);


  new TCanvas;

  // Let's repeat for the Icarus-only contour
  calc.SetL(kBaselineIcarus);

  Surface surfIcarus(&exptIcarus, &calc, kAxSinSq2ThetaMuMu, kAxDmSq);

  surfSBND.DrawContour(critSurf, kSolid, kRed);
  surfIcarus.DrawContour(critSurf, kSolid, kBlue);

  // MultiExperimentSBN sums the chisqs from its constituent parts. It also
  // makes sure to set the right baselines at the right times during the fit.
  MultiExperimentSBN exptMulti({&exptSBND, &exptIcarus}, {kSBND, kICARUS});
  Surface surfMulti(&exptMulti, &calc, kAxSinSq2ThetaMuMu, kAxDmSq);

  surfMulti.DrawContour(critSurf, kSolid, kMagenta);

  // Exercise - plot the MicroBooNE sensitivity and the three-detector joint
  // sensitivity.
}
