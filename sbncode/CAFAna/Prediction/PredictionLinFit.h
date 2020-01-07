#pragma once

#include "CAFAna/Prediction/IPrediction.h"

#include "TVectorD.h"

namespace ana
{
  /// Parameterize a collection of universes as a linear function of the syst
  /// knobs. Profiling over this should produce the same results as a
  /// covariance matrix fit.
  class PredictionLinFit: public IPrediction
  {
  public:
    PredictionLinFit(const std::vector<const ISyst*>& systs,
                     const IPrediction* pnom,
                     const std::vector<std::pair<SystShifts, const IPrediction*>>& univs);

    ~PredictionLinFit();

    Spectrum Predict(osc::IOscCalculator* calc) const override;
    Spectrum PredictSyst(osc::IOscCalculator* calc,
                         const SystShifts& syst) const override;

    Spectrum PredictComponent(osc::IOscCalculator* calc,
                              Flavors::Flavors_t flav,
                              Current::Current_t curr,
                              Sign::Sign_t sign) const override;
    Spectrum PredictComponentSyst(osc::IOscCalculator* calc,
                                  const SystShifts& syst,
                                  Flavors::Flavors_t flav,
                                  Current::Current_t curr,
                                  Sign::Sign_t sign) const override;

    void SaveTo(TDirectory* dir) const override;
    static std::unique_ptr<PredictionLinFit> LoadFrom(TDirectory* dir);

    protected:
      Ratio GetRatio(const SystShifts& shift) const;

      const std::vector<const ISyst*> fSysts;
      const IPrediction* fNom;
      std::vector<std::pair<SystShifts, const IPrediction*>> fUnivs;

      std::vector<TVectorD> fCoeffs;
  };
}
