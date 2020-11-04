#pragma once

#include "CAFAna/Core/SystShifts.h"
#include "CAFAna/Core/LoadFromFile.h"

#include <unordered_map>

class TDirectory;

namespace osc
{
  template<class T> class _IOscCalc;
  typedef _IOscCalc<double> IOscCalc;
  template<class T> class _IOscCalcAdjustable;
  typedef _IOscCalcAdjustable<double> IOscCalcAdjustable;
}

namespace ana
{
  /// Base class defining interface for experiments
  class IExperiment
  {
  public:
    virtual ~IExperiment() {}
    virtual double ChiSq(osc::IOscCalcAdjustable* osc,
                         const SystShifts& syst = SystShifts::Nominal()) const = 0;

    virtual void SaveTo(TDirectory* dir) const;
  };
}
