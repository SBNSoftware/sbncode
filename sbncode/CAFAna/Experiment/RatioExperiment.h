#pragma once

#include "CAFAna/Core/Spectrum.h"

#include "CAFAna/Experiment/IExperiment.h"
#include "CAFAna/Analysis/ExpInfo.h"

#include "OscLib/IOscCalculator.h"

#include <memory>
#include <vector>

namespace ana
{
  class RatioExperiment: public IExperiment
  {
  public:
    RatioExperiment(const IPrediction* predND,
                    const IPrediction* predFD,
                    const Spectrum& spectND,
                    const Spectrum& spectFD)
      : fPredND(predND), fPredFD(predFD), fSpectND(spectND), fSpectFD(spectFD)
    {
    }

    virtual double ChiSq(osc::IOscCalculatorAdjustable* osc,
                         const SystShifts& syst = SystShifts::Nominal()) const override;

    //    virtual void SaveTo(TDirectory* dir) const override;
    //    static std::unique_ptr<RatioExperiment> LoadFrom(TDirectory* dir);

  protected:
    const IPrediction* fPredND;
    const IPrediction* fPredFD;
    Spectrum fSpectND;
    Spectrum fSpectFD;
  };
}
