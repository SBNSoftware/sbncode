#include "CAFAna/Prediction/PredictionSBNExtrap.h"

#include "CAFAna/Core/LoadFromFile.h"
#include "CAFAna/Core/Ratio.h"

#include "CAFAna/Analysis/ExpInfo.h"

#include "OscLib/func/IOscCalculator.h"

#include "TDirectory.h"
#include "TObjString.h"
#include "TH1D.h"

namespace ana
{
  //----------------------------------------------------------------------
  // TODO need data and MC shifts/weights?
  PredictionSBNExtrap::PredictionSBNExtrap(Loaders& loadersND,
                                           Loaders& loadersFD,
                                           const HistAxis& axis,
                                           const Cut& cut,
                                           const SystShifts& shift,
                                           const Var& wei)
    : fPredND(loadersND, axis, cut, shift, wei),
      fPredFD(loadersFD, axis, cut, shift, wei),
      // TODO how on earth do fake ND oscillations get in here?
      fDataND(loadersND.GetLoader(Loaders::kData), axis, cut, shift, wei)
  {
  }

  //----------------------------------------------------------------------
  PredictionSBNExtrap::~PredictionSBNExtrap()
  {
  }

  //----------------------------------------------------------------------
  Spectrum PredictionSBNExtrap::Predict(osc::IOscCalculator* calc) const
  {
    return PredictComponent(calc,
                            Flavors::kAll,
                            Current::kBoth,
                            Sign::kBoth);
  }

  //----------------------------------------------------------------------
  Spectrum PredictionSBNExtrap::PredictComponent(osc::IOscCalculator* calc,
                                                 Flavors::Flavors_t flav,
                                                 Current::Current_t curr,
                                                 Sign::Sign_t sign) const
  {
    // TODO think hard about what ratios exactly we want to take
    // TODO the ND POT never comes into this :(
    ((osc::IOscCalculatorAdjustable*)calc)->SetL(kBaselineSBND);
    const Ratio r = fDataND / fPredND.Predict(calc);
    ((osc::IOscCalculatorAdjustable*)calc)->SetL(kBaselineIcarus);
    return r * fPredFD.PredictComponent(calc, flav, curr, sign);
  }

  //----------------------------------------------------------------------
  OscillatableSpectrum PredictionSBNExtrap::ComponentCC(int from, int to) const
  {
    /*
    if(from == +12 && to == +12) return fSBNExtrap->NueSurvComponent();
    if(from == -12 && to == -12) return fSBNExtrap->AntiNueSurvComponent();

    if(from == +12 && to == +14) return fSBNExtrap->NumuAppComponent();
    if(from == -12 && to == -14) return fSBNExtrap->AntiNumuAppComponent();

    if(from == +12 && to == +16) return fSBNExtrap->TauFromEComponent();
    if(from == -12 && to == -16) return fSBNExtrap->AntiTauFromEComponent();

    if(from == +14 && to == +12) return fSBNExtrap->NueAppComponent();
    if(from == -14 && to == -12) return fSBNExtrap->AntiNueAppComponent();

    if(from == +14 && to == +14) return fSBNExtrap->NumuSurvComponent();
    if(from == -14 && to == -14) return fSBNExtrap->AntiNumuSurvComponent();

    if(from == +14 && to == +16) return fSBNExtrap->TauFromMuComponent();
    if(from == -14 && to == -16) return fSBNExtrap->AntiTauFromMuComponent();
    */
    assert(0 && "Not reached");
  }

  //----------------------------------------------------------------------
  // Spectrum PredictionSBNExtrap::ComponentNC() const
  // {
  //   return fSBNExtrap->NCComponent();
  // }

  //----------------------------------------------------------------------
  void PredictionSBNExtrap::SaveTo(TDirectory* dir) const
  {
    TDirectory* tmp = gDirectory;

    dir->cd();

    TObjString("PredictionSBNExtrap").Write("type");

    fPredND.SaveTo(dir->mkdir("predND"));
    fPredFD.SaveTo(dir->mkdir("predFD"));
    fDataND.SaveTo(dir->mkdir("dataND"));

    tmp->cd();
  }

  //----------------------------------------------------------------------
  std::unique_ptr<PredictionSBNExtrap> PredictionSBNExtrap::LoadFrom(TDirectory* dir)
  {
    assert(dir->GetDirectory("predND"));
    assert(dir->GetDirectory("predFD"));
    assert(dir->GetDirectory("dataND"));

    // TODO are these leaks?
    return std::unique_ptr<PredictionSBNExtrap>(new PredictionSBNExtrap(
      *ana::LoadFrom<PredictionNoExtrap>(dir->GetDirectory("predND")).release(),
      *ana::LoadFrom<PredictionNoExtrap>(dir->GetDirectory("predFD")).release(),
      *ana::LoadFrom<Spectrum>(dir->GetDirectory("dataND")).release()));
  }
}
