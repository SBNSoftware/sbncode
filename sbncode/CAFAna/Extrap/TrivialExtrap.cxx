#include "CAFAna/Extrap/TrivialExtrap.h"

#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Core/Loaders.h"
#include "CAFAna/Core/SpectrumLoader.h"

#include "TDirectory.h"
#include "TObjString.h"

namespace ana
{
  //----------------------------------------------------------------------
  TrivialExtrap::TrivialExtrap(SpectrumLoaderBase& loaderNonswap,
                               SpectrumLoaderBase& loaderNue,
                               SpectrumLoaderBase& loaderNuTau,
                               SpectrumLoaderBase& loaderIntrinsic,
                               const HistAxis& axis,
                               const Cut& cut,
                               const SystShifts& shift,
                               const Var& wei)
    :
    fNueApp       (loaderNue,     axis, cut && kIsNueApp    && !kIsAntiNu, shift, wei),
    fNueAppAnti   (loaderNue,     axis, cut && kIsNueApp    &&  kIsAntiNu, shift, wei),

    fNumuSurv     (loaderNonswap, axis, cut && kIsNumuCC    && !kIsAntiNu, shift, wei),
    fNumuSurvAnti (loaderNonswap, axis, cut && kIsNumuCC    &&  kIsAntiNu, shift, wei),

    fNumuApp      (loaderNuTau,   axis, cut && kIsNumuApp   && !kIsAntiNu, shift, wei),
    fNumuAppAnti  (loaderNuTau,   axis, cut && kIsNumuApp   &&  kIsAntiNu, shift, wei),

    fNueSurv      (loaderIntrinsic, axis, cut && kIsBeamNue   && !kIsAntiNu, shift, wei),
    fNueSurvAnti  (loaderIntrinsic, axis, cut && kIsBeamNue   &&  kIsAntiNu, shift, wei),

    fTauFromE     (loaderNue,     axis, cut && kIsTauFromE  && !kIsAntiNu, shift, wei),
    fTauFromEAnti (loaderNue,     axis, cut && kIsTauFromE  &&  kIsAntiNu, shift, wei),

    fTauFromMu    (loaderNuTau,   axis, cut && kIsTauFromMu && !kIsAntiNu, shift, wei),
    fTauFromMuAnti(loaderNuTau,   axis, cut && kIsTauFromMu &&  kIsAntiNu, shift, wei),

    fNCFromNumu   (loaderNonswap, axis, cut && kIsNCFromNumu,     shift, wei),
    fNCFromNue    (loaderNonswap, axis, cut && kIsNCFromNue,      shift, wei)
  {
    // All swapped files are equally valid as a source of NCs. This
    // approximately doubles/triples our statistics. SpectrumLoader just adds
    // events and POT for both cases, which is the right thing to do.

    loaderNue  .AddReweightableSpectrum(fNCFromNumu, axis.GetVar1D(), kTrueLOverE, cut && kIsNCFromNumu, shift, wei);
    loaderNuTau.AddReweightableSpectrum(fNCFromNumu, axis.GetVar1D(), kTrueLOverE, cut && kIsNCFromNumu, shift, wei);

    loaderNue  .AddReweightableSpectrum(fNCFromNue, axis.GetVar1D(), kTrueLOverE, cut && kIsNCFromNue, shift, wei);
    loaderNuTau.AddReweightableSpectrum(fNCFromNue, axis.GetVar1D(), kTrueLOverE, cut && kIsNCFromNue, shift, wei);

    //Also load in intrinsic nues from nonswap file
    loaderNonswap.AddReweightableSpectrum(fNueSurv, axis.GetVar1D(), kTrueLOverE, cut && kIsBeamNue && !kIsAntiNu, shift, wei);
    loaderNonswap.AddReweightableSpectrum(fNueSurvAnti, axis.GetVar1D(), kTrueLOverE, cut && kIsBeamNue && kIsAntiNu, shift, wei);

  }


  //----------------------------------------------------------------------
  TrivialExtrap::TrivialExtrap(SpectrumLoaderBase& loaderNonswap,
                               SpectrumLoaderBase& loaderNue,
                               SpectrumLoaderBase& loaderNuTau,
                               SpectrumLoaderBase& loaderIntrinsic,
                               std::string label,
                               const Binning& bins,
                               const Var& var,
                               const Cut& cut,
                               const SystShifts& shift,
                               const Var& wei)
    :
    TrivialExtrap(loaderNonswap, loaderNue, loaderNuTau, loaderIntrinsic,
                  HistAxis(label, bins, var),
                  cut, shift, wei)
  {
  }

  //----------------------------------------------------------------------
  TrivialExtrap::TrivialExtrap(Loaders& loaders,
                               std::string label,
                               const Binning& bins,
                               const Var& var,
                               const Cut& cut,
                               const SystShifts& shift,
                               const Var& wei)
    : TrivialExtrap(loaders, HistAxis(label, bins, var), cut, shift, wei)
  {
  }

  //----------------------------------------------------------------------
  TrivialExtrap::TrivialExtrap(Loaders& loaders,
                               const HistAxis& axis,
                               const Cut& cut,
                               const SystShifts& shift,
                               const Var& wei)
    : TrivialExtrap(loaders.GetLoader(Loaders::kMC, Loaders::kNonSwap),
                    loaders.GetLoader(Loaders::kMC, Loaders::kNueSwap),
                    loaders.GetLoader(Loaders::kMC, Loaders::kNuTauSwap),
		    loaders.GetLoader(Loaders::kMC, Loaders::kIntrinsic),
                    axis, cut, shift, wei)
  {
  }

  //----------------------------------------------------------------------
  void TrivialExtrap::SaveTo(TDirectory* dir) const
  {
    TDirectory* tmp = gDirectory;

    dir->cd();

    TObjString("TrivialExtrap").Write("type");

    fNueApp.SaveTo(dir, "nue_app");
    fNueAppAnti.SaveTo(dir, "nue_app_anti");
    fNCFromNumu.SaveTo(dir, "nc_from_numu");
    fNCFromNue.SaveTo(dir, "nc_from_nue");
    fNumuSurv.SaveTo(dir, "numu_surv");
    fNumuSurvAnti.SaveTo(dir, "numu_surv_anti");
    fNumuApp.SaveTo(dir, "numu_app");
    fNumuAppAnti.SaveTo(dir, "numu_app_anti");
    fNueSurv.SaveTo(dir, "nue_surv");
    fNueSurvAnti.SaveTo(dir, "nue_surv_anti");
    fTauFromE.SaveTo(dir, "nutau_from_nue");
    fTauFromEAnti.SaveTo(dir, "nutau_from_nue_anti");
    fTauFromMu.SaveTo(dir, "nutau_from_numu");
    fTauFromMuAnti.SaveTo(dir, "nutau_from_numu_anti");

    tmp->cd();
  }

  //----------------------------------------------------------------------
  std::unique_ptr<TrivialExtrap> TrivialExtrap::LoadFrom(TDirectory* dir)
  {
    std::unique_ptr<TrivialExtrap> ret(new TrivialExtrap);

    // This is a lot of repetitive typing. Define some macros
#define LOAD_OSC(FIELD, LABEL) ret->FIELD = *OscillatableSpectrum::LoadFrom(dir, LABEL);

    ret->fNueApp        = *OscillatableSpectrum::LoadFrom(dir, "nue_app");
    ret->fNueAppAnti    = *OscillatableSpectrum::LoadFrom(dir, "nue_app_anti");
    ret->fNumuSurv      = *OscillatableSpectrum::LoadFrom(dir, "numu_surv");
    ret->fNumuSurvAnti  = *OscillatableSpectrum::LoadFrom(dir, "numu_surv_anti");
    ret->fNumuApp       = *OscillatableSpectrum::LoadFrom(dir, "numu_app");
    ret->fNumuAppAnti   = *OscillatableSpectrum::LoadFrom(dir, "numu_app_anti");
    ret->fNueSurv       = *OscillatableSpectrum::LoadFrom(dir, "nue_surv");
    ret->fNueSurvAnti   = *OscillatableSpectrum::LoadFrom(dir, "nue_surv_anti");
    ret->fTauFromE      = *OscillatableSpectrum::LoadFrom(dir, "nutau_from_nue");
    ret->fTauFromEAnti  = *OscillatableSpectrum::LoadFrom(dir, "nutau_from_nue_anti");
    ret->fTauFromMu     = *OscillatableSpectrum::LoadFrom(dir, "nutau_from_numu");
    ret->fTauFromMuAnti = *OscillatableSpectrum::LoadFrom(dir, "nutau_from_numu_anti");
    ret->fNCFromNumu    = *OscillatableSpectrum::LoadFrom(dir, "nc_from_numu");
    ret->fNCFromNue     = *OscillatableSpectrum::LoadFrom(dir, "nc_from_nue");

    return ret;
  }
}
