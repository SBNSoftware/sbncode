#if 0 // temporarily disable

#pragma once

#include "CAFAna/Prediction/IPrediction.h"

namespace ana
{
  class Loaders;
  class IPredictionGenerator;

  /// Parameterize a collection of universes as a function of the syst knobs
  class PredictionLinFit: public IPrediction
  {
  public:
    /// Direct creation from an ensemble of universes
    PredictionLinFit(const std::vector<const ISyst*>& systs,
                     const IPrediction* pnom,
                     const std::vector<std::pair<SystShifts, const IPrediction*>>& univs);

    /// Constructor in the PredictionInterp style
    PredictionLinFit(const std::vector<const ISyst*>& systs,
                     const IPredictionGenerator& predGen,
                     Loaders& loaders,
                     int nUniv);

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

    void DebugPlot(const ISyst* syst,
                   osc::IOscCalculator* calc,
                   Flavors::Flavors_t flav = Flavors::kAll,
                   Current::Current_t curr = Current::kBoth,
                   Sign::Sign_t sign = Sign::kBoth) const;

    // If \a savePattern is not empty, print each pad. If it contains "%s" then
    // multiple files will be written, one per systematic.
    void DebugPlots(osc::IOscCalculator* calc,
		    const std::string& savePattern = "",
		    Flavors::Flavors_t flav = Flavors::kAll,
		    Current::Current_t curr = Current::kBoth,
		    Sign::Sign_t sign = Sign::kBoth) const;

    void DebugPlotColz(const ISyst* syst,
                       osc::IOscCalculator* calc,
                       Flavors::Flavors_t flav = Flavors::kAll,
                       Current::Current_t curr = Current::kBoth,
                       Sign::Sign_t sign = Sign::kBoth) const;

    void DebugPlotsColz(osc::IOscCalculator* calc,
                        const std::string& savePattern = "",
                        Flavors::Flavors_t flav = Flavors::kAll,
                        Current::Current_t curr = Current::kBoth,
                        Sign::Sign_t sign = Sign::kBoth) const;

    void SaveTo(TDirectory* dir) const override;
    static std::unique_ptr<PredictionLinFit> LoadFrom(TDirectory* dir);

  protected:
    void InitFits() const;

    /// Helper for InitFits()
    std::vector<double> InitFitsBin(const std::vector<std::vector<double>>& M,
                                    const std::vector<double>& ds,
                                    const std::vector<std::vector<double>>& coords) const;

    std::vector<double> GetCoords(const SystShifts& shift) const;

    Ratio GetRatio(const SystShifts& shift) const;

    const std::vector<const ISyst*> fSysts;
    const IPrediction* fNom;
    std::vector<std::pair<SystShifts, const IPrediction*>> fUnivs;

    mutable std::vector<std::vector<double>> fCoeffs;
  };
}

#endif
