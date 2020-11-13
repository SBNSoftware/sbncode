#pragma once

#include "CAFAna/Experiment/IExperiment.h"

#include <memory>

namespace ana
{
  /// Very simple model allowing inclusion of reactor constraints
  class ReactorExperiment : public IExperiment
  {
  public:
    ReactorExperiment(double bestFit, double sigma)
      : fBestFit(bestFit), fSigma(sigma)
    {
    }

    virtual double ChiSq(osc::IOscCalculatorAdjustable* osc,
                         const SystShifts& shift = SystShifts::Nominal()) const override;

    void SaveTo(TDirectory* dir) const override;
    static std::unique_ptr<ReactorExperiment> LoadFrom(TDirectory* dir);
  protected:
    double fBestFit, fSigma;
  };

  /// A \ref ReactorExperiment initialized with the Nu2014 Daya Bay constraints
  const ReactorExperiment* DayaBayConstraint2014();

  /// Weighted average of all experiments as of first nue paper writing
  const ReactorExperiment* WorldReactorConstraint2015();

  /// Updated value for SecondAna based on the latest PDG
  const ReactorExperiment* WorldReactorConstraint2016();

  /// Reactor constraint from PDG 2017 update 
  const ReactorExperiment* WorldReactorConstraint2017();
}
