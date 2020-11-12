#pragma once

namespace osc{class IOscCalculatorAdjustable;}

#include "CAFAna/Core/SystShifts.h"
#include "CAFAna/Core/LoadFromFile.h"

#include <unordered_map>

class TDirectory;

namespace ana
{
  /// Base class defining interface for experiments
  class IExperiment
  {
  public:
    virtual ~IExperiment() {}
    virtual double ChiSq(osc::IOscCalculatorAdjustable* osc,
                         const SystShifts& syst = SystShifts::Nominal()) const = 0;

    virtual void SaveTo(TDirectory* dir) const;
  };
}
