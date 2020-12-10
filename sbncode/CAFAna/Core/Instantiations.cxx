#include "CAFAna/Core/SpectrumLoaderBase.h"

#include "CAFAna/Core/Spectrum.h"

//#include "CAFAna/Core/SpectrumConstructors.txx"

namespace ana
{
  // This is not what SpectrumConstructors.txx has, but this is what we need
  template<> Spectrum::Spectrum(SpectrumLoaderBase& loader,
                                const HistAxis& axis,
                                const Cut& cut,
                                const SystShifts& shift,
                                const Var& wei,
                                Spectrum::ESparse sparse)
    : Spectrum(LabelsAndBins(axis.GetLabels(), axis.GetBinnings()), sparse)
  {
    if(axis.HasVars()) loader.AddSpectrum(*this, axis.GetVar1D(), kNoSpillCut, cut, shift, wei);
  }

  // TODO why can't these be compiled into cafanacore now?
  //----------------------------------------------------------------------
  Spectrum::Spectrum(Spectrum&& rhs):
    fHist(std::move(rhs.fHist)),
    fPOT(rhs.fPOT),
    fLivetime(rhs.fLivetime),
    fAxis(rhs.fAxis)
  {
    std::swap(fReferences, rhs.fReferences);
    for(Spectrum** ref: fReferences) *ref = this;
  }

  //----------------------------------------------------------------------
  Spectrum& Spectrum::operator=(Spectrum&& rhs)
  {
    if(this == &rhs) return *this;

    fHist = std::move(rhs.fHist);
    fPOT = rhs.fPOT;
    fLivetime = rhs.fLivetime;
    fAxis = rhs.fAxis;

    std::swap(fReferences, rhs.fReferences);
    for(Spectrum** ref: fReferences) *ref = this;

    return *this;
  }

  /*
  template Spectrum::Spectrum(const std::string& label,
                              const Binning& bins,
                              SpectrumLoaderBase& loader,
                              const Var& var,
                              const Cut& cut,
                              const SystShifts& shift,
                              const Var& wei,
                              Spectrum::ESparse sparse);

  template Spectrum::Spectrum(SpectrumLoaderBase& loader,
                              const HistAxis& xAxis,
                              const HistAxis& yAxis,
                              const Cut& cut,
                              const SystShifts& shift,
                              const Var& wei,
                              ESparse sparse);

  template Spectrum::Spectrum(const std::string& xLabel,
                              const std::string& yLabel,
                              SpectrumLoaderBase& loader,
                              const Binning& binsx, const Var& varx,
                              const Binning& binsy, const Var& vary,
                              const Cut& cut,
                              const SystShifts& shift,
                              const Var& wei,
                              ESparse sparse);

  template Spectrum::Spectrum(SpectrumLoaderBase& loader,
                              const HistAxis& xAxis,
                              const HistAxis& yAxis,
                              const HistAxis& zAxis,
                              const Cut& cut,
                              const SystShifts& shift,
                              const Var& wei,
                              ESparse sparse);

  template Spectrum::Spectrum(const std::string& xLabel,
                              const std::string& yLabel,
                              const std::string& zLabel,
                              SpectrumLoaderBase& loader,
                              const Binning& binsx, const Var& varx,
                              const Binning& binsy, const Var& vary,
                              const Binning& binsz, const Var& varz,
                              const Cut& cut,
                              const SystShifts& shift,
                              const Var& wei,
                              ESparse sparse);

  template Spectrum::Spectrum(SpectrumLoaderBase& loader,
                              const _HistAxis<MultiVar>& axis,
                              const Cut& cut,
                              const SystShifts& shift,
                              const Var& wei);
  */
}

#include "CAFAna/Core/ReweightableSpectrum.h"
#include "CAFAna/Core/ReweightableSpectrumConstructors.txx"

namespace ana
{
  template ReweightableSpectrum::ReweightableSpectrum(SpectrumLoaderBase& loader,
                                                      const HistAxis& recoAxis,
                                                      const HistAxis& trueAxis,
                                                      const Cut& cut,
                                                      const SystShifts& shift,
                                                      const Var& wei);
}
