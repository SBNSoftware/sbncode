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
                               const HistAxis& axis,
                               const Cut& cut,
                               const SystShifts& shift,
                               const Var& wei)
    :
    fNueApp       (loaderNue,     axis, cut && kIsSig       && !kIsAntiNu, shift, wei),
    fNueAppAnti   (loaderNue,     axis, cut && kIsSig       &&  kIsAntiNu, shift, wei),

    fNumuSurv     (loaderNonswap, axis, cut && kIsNumuCC    && !kIsAntiNu, shift, wei),
    fNumuSurvAnti (loaderNonswap, axis, cut && kIsNumuCC    &&  kIsAntiNu, shift, wei),

    fNumuApp      (loaderNuTau,   axis, cut && kIsNumuApp   && !kIsAntiNu, shift, wei),
    fNumuAppAnti  (loaderNuTau,   axis, cut && kIsNumuApp   &&  kIsAntiNu, shift, wei),

    fNueSurv      (loaderNonswap, axis, cut && kIsBeamNue   && !kIsAntiNu, shift, wei),
    fNueSurvAnti  (loaderNonswap, axis, cut && kIsBeamNue   &&  kIsAntiNu, shift, wei),

    fTauFromE     (loaderNue,     axis, cut && kIsTauFromE  && !kIsAntiNu, shift, wei),
    fTauFromEAnti (loaderNue,     axis, cut && kIsTauFromE  &&  kIsAntiNu, shift, wei),

    fTauFromMu    (loaderNuTau,   axis, cut && kIsTauFromMu && !kIsAntiNu, shift, wei),
    fTauFromMuAnti(loaderNuTau,   axis, cut && kIsTauFromMu &&  kIsAntiNu, shift, wei),

    fNCFromNumu   (loaderNonswap, axis, cut && kIsNC && kIsNCFromNumu,     shift, wei),
    fNCFromNue    (loaderNonswap, axis, cut && kIsNC && kIsNCFromNue,      shift, wei)
  {
    // All swapped files are equally valid as a source of NCs. This
    // approximately doubles/triples our statistics. SpectrumLoader just adds
    // events and POT for both cases, which is the right thing to do.

    loaderNue  .AddReweightableSpectrum(fNCFromNumu, axis.GetMultiDVar(), cut && kIsNC && kIsNCFromNumu, shift, wei);
    loaderNuTau.AddReweightableSpectrum(fNCFromNumu, axis.GetMultiDVar(), cut && kIsNC && kIsNCFromNumu, shift, wei);

    loaderNue  .AddReweightableSpectrum(fNCFromNue, axis.GetMultiDVar(), cut && kIsNC && kIsNCFromNue, shift, wei);
    loaderNuTau.AddReweightableSpectrum(fNCFromNue, axis.GetMultiDVar(), cut && kIsNC && kIsNCFromNue, shift, wei);

  }


  //----------------------------------------------------------------------
  TrivialExtrap::TrivialExtrap(SpectrumLoaderBase& loaderNonswap,
                               SpectrumLoaderBase& loaderNue,
                               SpectrumLoaderBase& loaderNuTau,
                               std::string label,
                               const Binning& bins,
                               const Var& var,
                               const Cut& cut,
                               const SystShifts& shift,
                               const Var& wei)
    :
    TrivialExtrap(loaderNonswap, loaderNue, loaderNuTau,
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
    : TrivialExtrap(loaders.GetLoader(Loaders::kMC, ana::kBeam, Loaders::kNonSwap),
                    loaders.GetLoader(Loaders::kMC, ana::kBeam, Loaders::kNueSwap),
                    loaders.GetLoader(Loaders::kMC, ana::kBeam, Loaders::kNuTauSwap),
                    axis, cut, shift, wei)
  {
  }

  //----------------------------------------------------------------------
  void TrivialExtrap::SaveTo(TDirectory* dir) const
  {
    TDirectory* tmp = gDirectory;

    dir->cd();

    TObjString("TrivialExtrap").Write("type");

    fNueApp.SaveTo(dir->mkdir("nue_app"));
    fNueAppAnti.SaveTo(dir->mkdir("nue_app_anti"));
    fNCFromNumu.SaveTo(dir->mkdir("nc_from_numu"));
    fNCFromNue.SaveTo(dir->mkdir("nc_from_nue"));
    fNumuSurv.SaveTo(dir->mkdir("numu_surv"));
    fNumuSurvAnti.SaveTo(dir->mkdir("numu_surv_anti"));
    fNumuApp.SaveTo(dir->mkdir("numu_app"));
    fNumuAppAnti.SaveTo(dir->mkdir("numu_app_anti"));
    fNueSurv.SaveTo(dir->mkdir("nue_surv"));
    fNueSurvAnti.SaveTo(dir->mkdir("nue_surv_anti"));
    fTauFromE.SaveTo(dir->mkdir("nutau_from_nue"));
    fTauFromEAnti.SaveTo(dir->mkdir("nutau_from_nue_anti"));
    fTauFromMu.SaveTo(dir->mkdir("nutau_from_numu"));
    fTauFromMuAnti.SaveTo(dir->mkdir("nutau_from_numu_anti"));

    tmp->cd();
  }

  //----------------------------------------------------------------------
  std::unique_ptr<TrivialExtrap> TrivialExtrap::LoadFrom(TDirectory* dir)
  {
    std::unique_ptr<TrivialExtrap> ret(new TrivialExtrap);

    // This is a lot of repetitive typing. Define some macros
#define LOAD_OSC(FIELD, LABEL) assert(dir->GetDirectory(LABEL)); ret->FIELD = *OscillatableSpectrum::LoadFrom(dir->GetDirectory(LABEL));
#define LOAD_SPECT(FIELD, LABEL) assert(dir->GetDirectory(LABEL)); ret->FIELD = *Spectrum::LoadFrom(dir->GetDirectory(LABEL));

    LOAD_OSC(fNueApp,        "nue_app");
    LOAD_OSC(fNueAppAnti,    "nue_app_anti");
    LOAD_OSC(fNumuSurv,      "numu_surv");
    LOAD_OSC(fNumuSurvAnti,  "numu_surv_anti");
    LOAD_OSC(fNumuApp,       "numu_app");
    LOAD_OSC(fNumuAppAnti,   "numu_app_anti");
    LOAD_OSC(fNueSurv,       "nue_surv");
    LOAD_OSC(fNueSurvAnti,   "nue_surv_anti");
    LOAD_OSC(fTauFromE,      "nutau_from_nue");
    LOAD_OSC(fTauFromEAnti,  "nutau_from_nue_anti");
    LOAD_OSC(fTauFromMu,     "nutau_from_numu");
    LOAD_OSC(fTauFromMuAnti, "nutau_from_numu_anti");
    LOAD_OSC(fNCFromNumu,    "nc_from_numu");
    LOAD_OSC(fNCFromNue,     "nc_from_nue");

    return ret;
  }
}
