#pragma once

#include "CAFAna/Prediction/PredictionNoExtrap.h"

//#include "CAFAna/Prediction/PredictionGenerator.h"

namespace ana
{
  class Loaders;

  /// Prediction summing detector and dirt components
  class PredictionIncDirt: public IPrediction
  {
  public:
    PredictionIncDirt(SpectrumLoaderBase& loaderNonswap,
                      SpectrumLoaderBase& loaderNue,
                      SpectrumLoaderBase& loaderNuTau,
                      SpectrumLoaderBase& loaderIntrinsic,
                      SpectrumLoaderBase& loaderDirt,
                      const HistAxis& axis,
                      const Cut& cut,
                      const SystShifts& shift = kNoShift,
                      const Var& wei = kUnweighted);

    PredictionIncDirt(Loaders& loaders,
                      SpectrumLoaderBase& loaderDirt,
                      const HistAxis& axis,
                      const Cut& cut,
                      const SystShifts& shift = kNoShift,
                      const Var& wei = kUnweighted);

    virtual ~PredictionIncDirt();

    static std::unique_ptr<PredictionIncDirt> LoadFrom(TDirectory* dir);
    virtual void SaveTo(TDirectory* dir) const override;

    Spectrum PredictDet(osc::IOscCalc* calc) const
    {
      return fDet.Predict(calc);
    }

    Spectrum PredictDirt(osc::IOscCalc* calc) const
    {
      return fDirt.Predict(calc);
    }

    Spectrum PredictComponentDet(osc::IOscCalc* calc,
                                 Flavors::Flavors_t flav,
                                 Current::Current_t curr,
                                 Sign::Sign_t sign) const
    {
      return fDet.PredictComponent(calc, flav, curr, sign);
    }

    Spectrum PredictComponentDirt(osc::IOscCalc* calc,
                                  Flavors::Flavors_t flav,
                                  Current::Current_t curr,
                                  Sign::Sign_t sign) const
    {
      return fDirt.PredictComponent(calc, flav, curr, sign);
    }
    
    virtual Spectrum Predict(osc::IOscCalc* calc) const override
    {
      return PredictDet(calc) + PredictDirt(calc);
    }

    virtual Spectrum PredictComponent(osc::IOscCalc* calc,
                                      Flavors::Flavors_t flav,
                                      Current::Current_t curr,
                                      Sign::Sign_t sign) const override
    {
      return (PredictComponentDet (calc, flav, curr, sign) +
              PredictComponentDirt(calc, flav, curr, sign));
    }

  protected:
    PredictionIncDirt(std::unique_ptr<PredictionNoExtrap>&& det,
                      std::unique_ptr<PredictionNoExtrap>&& dirt)
      : fDet(*det), fDirt(*dirt)
    {
    }

    PredictionNoExtrap fDet, fDirt;
  };

  // TODO how best to write a PredictionGenerator for this?
}
