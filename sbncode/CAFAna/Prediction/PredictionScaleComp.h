#pragma once

#include "CAFAna/Core/ISyst.h"
#include "CAFAna/Core/SystShifts.h"
#include "CAFAna/Prediction/IPrediction.h"
#include "CAFAna/Systs/SystComponentScale.h"

#include "OscLib/IOscCalc.h"

namespace ana
{
  /// \brief Prediction broken down into arbitrary components whose scales can
  /// be varied independently.
  class PredictionScaleComp : public IPrediction
  {
  public:
    /// \param cut Cut applied to all histograms
    /// \param truthcuts Prediction will be broken down into N components
    ///                  following these cuts.
    PredictionScaleComp(SpectrumLoaderBase& loader,
                        const HistAxis&     axis,
                        Cut                 cut,
                        const std::vector<const SystComponentScale*>& systs,
                        const SystShifts&   shift = kNoShift,
                        const Var&          wei = kUnweighted);
    /// Constructor to take two HistAxis's to weight 2D spectra
    PredictionScaleComp(SpectrumLoaderBase& loader,
                        const HistAxis&     axis1,
                        const HistAxis&     axis2,
                        Cut                 cut,
                        const std::vector<const SystComponentScale*>& systs,
                        const SystShifts&   shift = kNoShift,
                        const Var&          wei = kUnweighted);

    /// This is for the FD via PredictionNoExtrap
    PredictionScaleComp(SpectrumLoaderBase& loaderNonswap,
                        SpectrumLoaderBase& loaderNue,
                        SpectrumLoaderBase& loaderNuTau,
                        SpectrumLoaderBase& loaderIntrinsic,
                        const HistAxis&     axis,
                        Cut                 cut,
                        const std::vector<const SystComponentScale*>& systs,
                        const SystShifts&   shift = kNoShift,
                        const Var&          wei = kUnweighted);

    virtual ~PredictionScaleComp();

    virtual Spectrum Predict(osc::IOscCalc* osc) const override
    {
      return fTotal->Predict(osc);
    }

    virtual Spectrum PredictSyst(osc::IOscCalc* osc,
                                 const SystShifts&    syst) const override
    {
      return PredictComponentSyst(osc, syst,
                                  Flavors::kAll, Current::kBoth, Sign::kBoth);
    }

    virtual Spectrum PredictComponent(osc::IOscCalc* calc,
                                      Flavors::Flavors_t flav,
                                      Current::Current_t curr,
                                      Sign::Sign_t sign) const override
    {
      return fTotal->PredictComponent(calc, flav, curr, sign);
    }

    virtual Spectrum PredictComponentSyst(osc::IOscCalc* calc,
                                          const SystShifts& syst,
                                          Flavors::Flavors_t flav,
                                          Current::Current_t curr,
                                          Sign::Sign_t sign) const override;

    Spectrum PredictCategory(osc::IOscCalc* osc,
                             const SystComponentScale* syst) const;

    static std::unique_ptr<PredictionScaleComp> LoadFrom(TDirectory* dir);
    virtual void SaveTo(TDirectory* dir) const override;

  protected:
    PredictionScaleComp(const IPrediction* total,
                        const std::vector<const IPrediction*>& preds,
                        const std::vector<const SystComponentScale*>& systs);

    std::vector<const SystComponentScale*> fSysts;
    std::vector<const IPrediction*> fPreds;

    const IPrediction* fTotal;
  };
}
