#pragma once

#include "CAFAna/Core/FwdDeclare.h"
#include "CAFAna/Prediction/IPrediction.h"
#include "CAFAna/Experiment/IExperiment.h"

namespace ana
{
  /// Compare a single data spectrum to the MC expectation
  class SingleSampleExperiment: public IExperiment
  {
  public:
    /// \param pred   Source of oscillated MC beam predictions
    /// \param data   Data spectrum to compare to
    SingleSampleExperiment(const IPrediction* pred,
                           const Spectrum& data);

    virtual ~SingleSampleExperiment();

    virtual double ChiSq(osc::IOscCalcAdjustable* osc,
                         const SystShifts& syst = SystShifts::Nominal()) const override;

    /*
    stan::math::var LogLikelihood(osc::_IOscCalcAdjustable<stan::math::var> *osc,
                                  const SystShifts &syst = kNoShift) const override;
    */

    virtual void SaveTo(TDirectory* dir) const override;
    static std::unique_ptr<SingleSampleExperiment> LoadFrom(TDirectory* dir);

    // Didn't make provisions for copying fMC
    SingleSampleExperiment(const SingleSampleExperiment&) = delete;
    SingleSampleExperiment& operator=(const SingleSampleExperiment&) = delete;

    // need to explicitly declare move constructor since copy constructor is deleted
    SingleSampleExperiment(SingleSampleExperiment&& s)
      : fMC(s.fMC), fData(std::move(s.fData))
    {
      s.fMC = nullptr;
    }

    void SetMaskHist(double xmin=0, double xmax=-1, 
		     double ymin=0, double ymax=-1);

  protected:
    virtual void ApplyMask(Eigen::ArrayXd& a, Eigen::ArrayXd& b) const;

    const IPrediction* fMC;
    Spectrum fData;
    Eigen::ArrayXd fMaskA;
  };
}
