#pragma once

#include "CAFAna/Decomp/IDecomp.h"

#include "CAFAna/Core/Loaders.h"

namespace ana
{
  /// Uses MC for NC and \f$ \nu_e \f$ CC components, assigns remainder of data to \f$ \nu_\mu \f$ CC
  class NumuDecomp: public IDecomp
  {
  public:
    NumuDecomp(SpectrumLoaderBase& loaderMC,
               SpectrumLoaderBase& loaderData,
               const std::string& label,
               const Binning& bins,
               const Var& var,
               const Cut& cut,
               const SystShifts& shiftMC = kNoShift,
               const SystShifts& shiftData = kNoShift,
               const Var& wei = kUnweighted);

    NumuDecomp(SpectrumLoaderBase& loaderMC,
               SpectrumLoaderBase& loaderData,
               const HistAxis& axis,
               const Cut& cut,
               const SystShifts& shiftMC = kNoShift,
               const SystShifts& shiftData = kNoShift,
               const Var& wei = kUnweighted);

    NumuDecomp(Loaders& loaders,
               const HistAxis& axis,
               const Cut& cut,
               const SystShifts& shiftMC = kNoShift,
               const SystShifts& shiftData = kNoShift,
               const Var& wei = kUnweighted);

    Spectrum NumuComponent()     const override;
    Spectrum AntiNumuComponent() const override;
    Spectrum NCComponent()       const override {return fNC;}
    Spectrum NueComponent()      const override {return fNue;}
    Spectrum AntiNueComponent()  const override {return fAntiNue;}

    virtual Spectrum MC_NumuComponent() const {return fNumu;}
    virtual Spectrum MC_AntiNumuComponent() const {return fAntiNumu;}
    virtual Spectrum Data_Component () const {return fData;}

    void SaveTo(TDirectory* dir) const override;
    static std::unique_ptr<NumuDecomp> LoadFrom(TDirectory* dir);

  protected:
    NumuDecomp()
      : fData(0, {}, {}, 0, 0),
        fNC  (0, {}, {}, 0, 0),
        fNue (0, {}, {}, 0, 0),
        fAntiNue(0, {}, {}, 0, 0),
        fNumu(0, {}, {}, 0, 0),
        fAntiNumu(0, {}, {}, 0, 0),
        fNotNumu(0, {}, {}, 0, 0)
    {};

    Spectrum fData;

    Spectrum fNC;

    Spectrum fNue;
    Spectrum fAntiNue;

    Spectrum fNumu;
    Spectrum fAntiNumu;

    Spectrum fNotNumu;
  };
}
