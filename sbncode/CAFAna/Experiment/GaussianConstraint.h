#pragma once

#include "CAFAna/Experiment/IExperiment.h"
#include "CAFAna/Core/IFitVar.h"

#include <memory>

namespace ana
{
  /// A simple Gaussian constraint on an arbitrary IFitVar
  class GaussianConstraint : public IExperiment
  {
  public:
  GaussianConstraint(const IFitVar* var, double mean, double sigma)
    : fVar(var), fMean(mean), fSigma(sigma)
    {
    }

    virtual double ChiSq(osc::IOscCalcAdjustable* osc,
                         const SystShifts& shift = SystShifts::Nominal()) const override;

  protected:
    const IFitVar* fVar;
    double fMean, fSigma;
  };

}
