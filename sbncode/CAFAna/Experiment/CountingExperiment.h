#pragma once

#include "CAFAna/Experiment/IExperiment.h"
#include "CAFAna/Core/Spectrum.h"


namespace ana
{
  class IPrediction;

  /// Compare a data spectrum to MC expectation using only the event count
  class CountingExperiment: public IExperiment
  {
  public:
    CountingExperiment(const IPrediction* p, const Spectrum& d, const Spectrum& cosmic);
    /// Version without cosmics may be wanted for MC studies
    CountingExperiment(const IPrediction* p, const Spectrum& d) : fMC(p), fData(d), fCosmic(0) {}
    ~CountingExperiment();
    virtual double ChiSq(osc::IOscCalculatorAdjustable* osc,
                         const SystShifts& syst = SystShifts::Nominal()) const override;

    virtual void SaveTo(TDirectory* dir) const override;
    static std::unique_ptr<CountingExperiment> LoadFrom(TDirectory* dir);
  protected:
    const IPrediction* fMC;
    Spectrum fData;
    TH1* fCosmic;
  };
}
