#include "CAFAna/Decomp/NumuDecomp.h"

#include "CAFAna/Cuts/TruthCuts.h"
#include "CAFAna/Core/Ratio.h"

#include <cassert>

#include "TDirectory.h"
#include "TObjString.h"

namespace ana
{
  //----------------------------------------------------------------------
  NumuDecomp::NumuDecomp(SpectrumLoaderBase& loaderMC,
                         SpectrumLoaderBase& loaderData,
                         const std::string& label,
                         const Binning& bins,
                         const Var& var,
                         const Cut& cut,
                         const SystShifts& shiftMC,
                         const SystShifts& shiftData,
                         const Var& wei)
    : NumuDecomp(loaderMC, loaderData, HistAxis(label, bins, var),
                 cut, shiftMC, shiftData, wei)
  {
  }

  //----------------------------------------------------------------------
  NumuDecomp::NumuDecomp(SpectrumLoaderBase& loaderMC,
                         SpectrumLoaderBase& loaderData,
                         const HistAxis& axis,
                         const Cut& cut,
                         const SystShifts& shiftMC,
                         const SystShifts& shiftData,
                         const Var& wei)
    : fData    (loaderData, axis, cut,                         shiftData, wei),
      fNC      (loaderMC,   axis, cut && kIsNC,                  shiftMC, wei),
      fNue     (loaderMC,   axis, cut && kIsBeamNue&&!kIsAntiNu, shiftMC, wei),
      fAntiNue (loaderMC,   axis, cut && kIsBeamNue&& kIsAntiNu, shiftMC, wei),
      fNumu    (loaderMC,   axis, cut && kIsNumuCC &&!kIsAntiNu, shiftMC, wei),
      fAntiNumu(loaderMC,   axis, cut && kIsNumuCC && kIsAntiNu, shiftMC, wei),
      fNotNumu (loaderMC,   axis, cut && !kIsNumuCC,             shiftMC, wei)
  {
  }

  //----------------------------------------------------------------------
  NumuDecomp::NumuDecomp(Loaders& loaders,
                         const HistAxis& axis,
                         const Cut& cut,
                         const SystShifts& shiftMC,
                         const SystShifts& shiftData,
                         const Var& wei)
    : NumuDecomp(loaders.GetLoader(caf::kNEARDET, Loaders::kMC),
                 loaders.GetLoader(caf::kNEARDET, Loaders::kData),
                 axis, cut, shiftMC, shiftData, wei)
  {
  }

  //----------------------------------------------------------------------
  Spectrum NumuDecomp::NumuComponent() const
  {
    // Subtract backgrounds from data using MC, then split into nu/antinu
    // components using MC ratio
    return (fNumu/(fNumu+fAntiNumu))*(fData-fNotNumu);
  }

  //----------------------------------------------------------------------
  Spectrum NumuDecomp::AntiNumuComponent() const
  {
    return (fAntiNumu/(fNumu+fAntiNumu))*(fData-fNotNumu);
  }

  //----------------------------------------------------------------------
  void NumuDecomp::SaveTo(TDirectory* dir) const
  {
    TDirectory* tmp = gDirectory;

    dir->cd();

    TObjString("NumuDecomp").Write("type");

    fNC.SaveTo(dir->mkdir("nc_comp"));
    fData.SaveTo(dir->mkdir("data_comp"));
    fNue.SaveTo(dir->mkdir("nue_comp"));
    fAntiNue.SaveTo(dir->mkdir("antinue_comp"));
    fNumu.SaveTo(dir->mkdir("numu_comp"));
    fAntiNumu.SaveTo(dir->mkdir("antinumu_comp"));
    fNotNumu.SaveTo(dir->mkdir("notnumu_comp"));

    tmp->cd();
  }

  //----------------------------------------------------------------------
  std::unique_ptr<NumuDecomp> NumuDecomp::LoadFrom(TDirectory* dir)
  {
    std::unique_ptr<NumuDecomp> ret(new NumuDecomp);

    // This is a lot of repetitive typing. Define a macro
#define LOAD_SPECT(FIELD, LABEL) assert(dir->GetDirectory(LABEL)); ret->FIELD = *Spectrum::LoadFrom(dir->GetDirectory(LABEL));

    LOAD_SPECT(fNC,       "nc_comp");
    LOAD_SPECT(fData,     "data_comp");
    LOAD_SPECT(fNue,      "nue_comp");
    LOAD_SPECT(fAntiNue,  "antinue_comp");
    LOAD_SPECT(fNumu,     "numu_comp");
    LOAD_SPECT(fAntiNumu, "antinumu_comp");
    LOAD_SPECT(fNotNumu,  "notnumu_comp");

    return ret;
  }
}
