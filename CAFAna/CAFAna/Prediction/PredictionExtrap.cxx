#include "CAFAna/Prediction/PredictionExtrap.h"

#include "CAFAna/Extrap/IExtrap.h"
#include "CAFAna/Core/LoadFromFile.h"

#include "TDirectory.h"
#include "TObjString.h"
#include "TH1D.h"

namespace ana
{
  //----------------------------------------------------------------------
  PredictionExtrap::PredictionExtrap(IExtrap* extrap)
    : fExtrap(extrap)
  {
  }

  //----------------------------------------------------------------------
  PredictionExtrap::~PredictionExtrap()
  {
    //    delete fExtrap;
  }

  //----------------------------------------------------------------------
  Spectrum PredictionExtrap::Predict(osc::IOscCalculator* calc) const
  {
    return PredictComponent(calc,
                            Flavors::kAll,
                            Current::kBoth,
                            Sign::kBoth);
  }

  //----------------------------------------------------------------------
  Spectrum PredictionExtrap::PredictComponent(osc::IOscCalculator* calc,
                                              Flavors::Flavors_t flav,
                                              Current::Current_t curr,
                                              Sign::Sign_t sign) const
  {
    Spectrum ret = fExtrap->NCComponent(); // Get binning
    ret.Clear();

    if(curr & Current::kCC){
      if(flav & Flavors::kNuEToNuE    && sign & Sign::kNu)     {
	ret += fExtrap->NueSurvComponent().    Oscillated(calc, +12, +12);
      }
      if(flav & Flavors::kNuEToNuE    && sign & Sign::kAntiNu) ret += fExtrap->AntiNueSurvComponent().Oscillated(calc, -12, -12);

      if(flav & Flavors::kNuEToNuMu   && sign & Sign::kNu)     ret += fExtrap->NumuAppComponent().    Oscillated(calc, +12, +14);
      if(flav & Flavors::kNuEToNuMu   && sign & Sign::kAntiNu) ret += fExtrap->AntiNumuAppComponent().Oscillated(calc, -12, -14);

      if(flav & Flavors::kNuEToNuTau  && sign & Sign::kNu)     ret += fExtrap->TauFromEComponent().    Oscillated(calc, +12, +16);
      if(flav & Flavors::kNuEToNuTau  && sign & Sign::kAntiNu) ret += fExtrap->AntiTauFromEComponent().Oscillated(calc, -12, -16);

      if(flav & Flavors::kNuMuToNuE   && sign & Sign::kNu)     ret += fExtrap->NueAppComponent().    Oscillated(calc, +14, +12);
      if(flav & Flavors::kNuMuToNuE   && sign & Sign::kAntiNu) ret += fExtrap->AntiNueAppComponent().Oscillated(calc, -14, -12);

      if(flav & Flavors::kNuMuToNuMu  && sign & Sign::kNu)     ret += fExtrap->NumuSurvComponent().    Oscillated(calc, +14, +14);
      if(flav & Flavors::kNuMuToNuMu  && sign & Sign::kAntiNu) ret += fExtrap->AntiNumuSurvComponent().Oscillated(calc, -14, -14);

      if(flav & Flavors::kNuMuToNuTau && sign & Sign::kNu)     ret += fExtrap->TauFromMuComponent().    Oscillated(calc, +14, +16);
      if(flav & Flavors::kNuMuToNuTau && sign & Sign::kAntiNu) ret += fExtrap->AntiTauFromMuComponent().Oscillated(calc, -14, -16);
    }
    if(curr & Current::kNC){
      assert(flav == Flavors::kAll); // Don't know how to calculate anything else
      assert(sign == Sign::kBoth);   // Why would you want to split NCs out by sign?

      ret += fExtrap->NCComponent();
    }

    return ret;
  }

  //----------------------------------------------------------------------
  OscillatableSpectrum PredictionExtrap::ComponentCC(int from, int to) const
  {
    if(from == +12 && to == +12) return fExtrap->NueSurvComponent();
    if(from == -12 && to == -12) return fExtrap->AntiNueSurvComponent();

    if(from == +12 && to == +14) return fExtrap->NumuAppComponent();
    if(from == -12 && to == -14) return fExtrap->AntiNumuAppComponent();

    if(from == +12 && to == +16) return fExtrap->TauFromEComponent();
    if(from == -12 && to == -16) return fExtrap->AntiTauFromEComponent();

    if(from == +14 && to == +12) return fExtrap->NueAppComponent();
    if(from == -14 && to == -12) return fExtrap->AntiNueAppComponent();

    if(from == +14 && to == +14) return fExtrap->NumuSurvComponent();
    if(from == -14 && to == -14) return fExtrap->AntiNumuSurvComponent();

    if(from == +14 && to == +16) return fExtrap->TauFromMuComponent();
    if(from == -14 && to == -16) return fExtrap->AntiTauFromMuComponent();

    assert(0 && "Not reached");
  }

  //----------------------------------------------------------------------
  Spectrum PredictionExtrap::ComponentNC() const
  {
    return fExtrap->NCComponent();
  }

  //----------------------------------------------------------------------
  void PredictionExtrap::SaveTo(TDirectory* dir) const
  {
    TDirectory* tmp = gDirectory;

    dir->cd();

    TObjString("PredictionExtrap").Write("type");

    fExtrap->SaveTo(dir->mkdir("extrap"));

    tmp->cd();
  }

  //----------------------------------------------------------------------
  std::unique_ptr<PredictionExtrap> PredictionExtrap::LoadFrom(TDirectory* dir)
  {
    assert(dir->GetDirectory("extrap"));
    IExtrap* extrap = ana::LoadFrom<IExtrap>(dir->GetDirectory("extrap")).release();

    return std::unique_ptr<PredictionExtrap>(new PredictionExtrap(extrap));
  }
}
