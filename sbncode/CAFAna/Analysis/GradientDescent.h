#pragma once

#include "Minuit2/MnApplication.h"

namespace ana
{
  /// A minimalistic gradient descent fitter to complement MINUIT's more
  /// elaborate offerings
  class GradientDescent: public ROOT::Minuit2::MnApplication
  {
  public:
    GradientDescent(const ROOT::Minuit2::FCNGradientBase& func,
                    const ROOT::Minuit2::MnUserParameters& pars);

    virtual ROOT::Minuit2::FunctionMinimum operator()(unsigned int maxfcn,
                                                      double tolerance) override;

    virtual const ROOT::Minuit2::ModularFunctionMinimizer& Minimizer() const override
    {
      abort();
    }

  protected:
    ROOT::Minuit2::FunctionMinimum Package(const std::vector<double>& pt,
                                           double chi, int ncalls) const;

    double Magnitude(const std::vector<double>& xs) const;
    void MakeUnit(std::vector<double>& xs) const;

    const ROOT::Minuit2::FCNGradientBase& fFunc;
    const ROOT::Minuit2::MnUserParameters& fPars;
  };
}
