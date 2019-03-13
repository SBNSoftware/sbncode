#pragma once

#include "CAFAna/Experiment/IExperiment.h"
#include "CAFAna/Prediction/IPrediction.h"
#include "CAFAna/Core/Spectrum.h"

namespace ana
{
  /// Dummy syst to communicate with \ref SingleSampleExperiment
  class CosmicBkgScaleSyst: public ISyst
  {
  public:
    CosmicBkgScaleSyst() : ISyst("cosmicScale", "Cosmic background scale") {}
    void Shift(double, Restorer&, caf::StandardRecord*, double&) const {}
  };

  extern const CosmicBkgScaleSyst kCosmicBkgScaleSyst;

  /// Compare a single data spectrum to the MC + cosmics expectation
  class SingleSampleExperiment: public IExperiment
  {
  public:
    /// \param pred   Source of oscillated MC beam predictions
    /// \param data   Data spectrum to compare to
    /// \param cosmic Cosmic ray background component
    /// \param cosmicScaleError fractional uncertainty on cosmic normalization
    SingleSampleExperiment(const IPrediction* pred,
                           const Spectrum& data,
                           const Spectrum& cosmic,
                           double cosmicScaleError = 0);

    /// \brief Fallback to manual cosmic scaling
    ///
    /// \a cosmic must be already scaled so that its bin contents can be
    /// directly summed onto \a data. If you're using the out-of-time part of
    /// the beam spill, the easiest thing to do is to pass \ref
    /// kTimingSidebandWeight as the weight argument when you fill it.
    SingleSampleExperiment(const IPrediction* pred,
                           const Spectrum& data,
                           const TH1D* cosmic,
                           double cosmicScaleError = 0);

    /// In MC studies you might not want to bother with cosmics
    SingleSampleExperiment(const IPrediction* pred,
                           const Spectrum& data)
      : fMC(pred), fData(data), fCosmic(0), fMask(0)
    {
    }

    virtual ~SingleSampleExperiment();

    virtual double ChiSq(osc::IOscCalculatorAdjustable* osc,
                         const SystShifts& syst = SystShifts::Nominal()) const override;

    virtual void Derivative(osc::IOscCalculator* calc,
                            const SystShifts& shift,
                            std::unordered_map<const ISyst*, double>& dch) const override;

    virtual void SaveTo(TDirectory* dir) const override;
    static std::unique_ptr<SingleSampleExperiment> LoadFrom(TDirectory* dir);

    // Didn't make provisions for copying fCosmic or fMC
    SingleSampleExperiment(const SingleSampleExperiment&) = delete;
    SingleSampleExperiment& operator=(const SingleSampleExperiment&) = delete;

    // need to explicitly declare move constructor since copy constructor is deleted
    SingleSampleExperiment(SingleSampleExperiment&& s)
      : fMC(s.fMC), fData(std::move(s.fData)), fCosmic(s.fCosmic), fCosmicScaleError(s.fCosmicScaleError)
    {
      s.fMC = nullptr;
      s.fCosmic = nullptr;
      s.fCosmicScaleError = 0;
    };

    void SetMaskHist(double xmin=0, double xmax=-1, 
		     double ymin=0, double ymax=-1);

  protected:
    TH1D* PredHistIncCosmics(osc::IOscCalculator* calc,
                             const SystShifts& syst) const;

    const IPrediction* fMC;
    Spectrum fData;
    TH1D* fCosmic;
    TH1* fMask;
    
    double fCosmicScaleError;
  };
}
