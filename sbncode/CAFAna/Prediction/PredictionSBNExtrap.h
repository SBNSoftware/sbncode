#pragma once

#include "CAFAna/Prediction/IPrediction.h"

#include "CAFAna/Prediction/PredictionNoExtrap.h"

#include "CAFAna/Core/Loaders.h"

namespace ana
{
  /// TODO comment
  class PredictionSBNExtrap: public IPrediction
  {
  public:
    PredictionSBNExtrap(Loaders& loadersND,
                        Loaders& loadersFD,
                        const HistAxis& axis,
                        const Cut& cut,
                        const SystShifts& shift_mc = kNoShift,
                        const Var& wei_mc = kUnweighted,
                        const SystShifts& shift_data = kNoShift,
                        const Var& wei_data = kUnweighted);
    virtual ~PredictionSBNExtrap();

    virtual Spectrum Predict(osc::IOscCalculator* calc) const override;

    virtual Spectrum PredictComponent(osc::IOscCalculator* calc,
                                      Flavors::Flavors_t flav,
                                      Current::Current_t curr,
                                      Sign::Sign_t sign) const override;

    OscillatableSpectrum ComponentCC(int from, int to) const override;
    //    Spectrum ComponentNC() const override;

    virtual void SaveTo(TDirectory* dir) const override;
    static std::unique_ptr<PredictionSBNExtrap> LoadFrom(TDirectory* dir);

    PredictionSBNExtrap() = delete;

  protected:
    PredictionSBNExtrap(const PredictionNoExtrap& pn,
                        const PredictionNoExtrap& pf,
                        const Spectrum& d)
      : fPredND(pn), fPredFD(pf), fDataND(d)
    {
    }

    PredictionNoExtrap fPredND, fPredFD;
    Spectrum fDataND;
  };

  class SBNExtrapGenerator: public IPredictionGenerator
  {
  public:
    SBNExtrapGenerator(Loaders& loaders_nd,
                       const HistAxis& ax,
                       const Cut& cut,
                       const Var& wei_mc,
                       const SystShifts& shift_data = kNoShift,
                       const Var& wei_data = kUnweighted);

    std::unique_ptr<IPrediction> Generate(Loaders& loaders_fd,
                                          const SystShifts& shiftMC = kNoShift) const override;

  protected:
    Loaders& fLoadersND;
    HistAxis fAxis;
    Cut fCut;
    Var fWeightMC;
    SystShifts fShiftData;
    Var fWeightData;
  };
}
