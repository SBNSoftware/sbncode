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

    virtual double ChiSq(osc::IOscCalculatorAdjustable* osc,
                         const SystShifts& shift = SystShifts::Nominal()) const override;

    virtual void Derivative(osc::IOscCalculator* calc,
                            const SystShifts& shift,
                            std::unordered_map<const ISyst*, double>& dchi) const override
    {
      // Empty implementation (rather than clearing the vector) means we have
      // zero derivative wrt systematics.
    }

  protected:
    const IFitVar* fVar;
    double fMean, fSigma;
  };

}
