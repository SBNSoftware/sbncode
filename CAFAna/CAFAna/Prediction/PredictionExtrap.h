#pragma once

#include "CAFAna/Prediction/IPrediction.h"

namespace ana
{
  class IExtrap;

  /// Take the output of an extrapolation and oscillate it as required
  class PredictionExtrap: public IPrediction
  {
  public:
    /// Takes ownership of \a extrap
    PredictionExtrap(IExtrap* extrap);
    virtual ~PredictionExtrap();

    virtual Spectrum Predict(osc::IOscCalculator* calc) const override;

    virtual Spectrum PredictComponent(osc::IOscCalculator* calc,
                                      Flavors::Flavors_t flav,
                                      Current::Current_t curr,
                                      Sign::Sign_t sign) const override;

    OscillatableSpectrum ComponentCC(int from, int to) const override;
    Spectrum ComponentNC() const override;

    virtual void SaveTo(TDirectory* dir) const override;
    static std::unique_ptr<PredictionExtrap> LoadFrom(TDirectory* dir);

    PredictionExtrap() = delete;

    IExtrap* GetExtrap() const {return fExtrap;}
  protected:

    IExtrap* fExtrap;
  };
}
