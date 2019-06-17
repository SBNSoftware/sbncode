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

    virtual void Derivative(osc::IOscCalculator* calc,
                            const SystShifts& shift,
                            std::unordered_map<const ISyst*, double>& dchi) const
    {
      // Optional to implement
      //
      // NB the convention is to *add* your contribution to the dchi values.
      //
      // If unimplemented, this default will be called, signaling no result to
      // the caller
      dchi.clear();
    }

    virtual void SaveTo(TDirectory* dir) const;
  };
}
