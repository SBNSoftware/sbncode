#pragma once

#include "CAFAna/Experiment/IExperiment.h"
#include "OscLib/func/IOscCalculator.h"

#include <memory>
#include <vector>

namespace ana
{
  /// Combine multiple component experiments
  class MultiExperimentSBN: public IExperiment
  {
  public:
    // We may want to rewrite this as a pair<extp, L>, but this works for now
    // also, to prevent mistakes it's probably best to define
    // kSBND, kUBoone and kIcarus somewhere. But again, works for now
    MultiExperimentSBN(std::vector<const IExperiment*> expts = {},
                       std::vector<float> baselines = {}) 
        : fExpts(expts), fLs(baselines)
    {
      fSystCorrelations.resize(expts.size());
    }

    void Add(const IExperiment* expt, float l){
                  fExpts.push_back(expt);
                  fLs.push_back(l);
                  }

    virtual double ChiSq(osc::IOscCalculatorAdjustable* osc,
                         const SystShifts& syst = SystShifts::Nominal()) const override;

    virtual void Derivative(osc::IOscCalculator* calc,
                            const SystShifts& shift,
                            std::unordered_map<const ISyst*, double>& dchi) const override;

    /// For the subexperiment \a idx, set up a mapping between systematics
    ///
    /// Each element in the vector is a pair from a "primary" systematic to a
    /// "secondary". When this MultiExperiment is called with a primary
    /// systematic shifted, the sub-experiment will be called with the
    /// secondary systematic set to the same value (and the primary unset).
    ///
    /// You can pass NULL for a secondary to indicate that the systematic
    /// simply has no effect on the experiment in question and should be
    /// filtered out.
    void SetSystCorrelations(int idx,
                             const std::vector<std::pair<const ISyst*,
                                                         const ISyst*>>& corrs);

    virtual void SaveTo(TDirectory* dir) const override;
    static std::unique_ptr<MultiExperimentSBN> LoadFrom(TDirectory* dir);

  protected:
    std::vector<const IExperiment*> fExpts;
    std::vector<float> fLs;
    std::vector<std::vector<std::pair<const ISyst*, const ISyst*>>> fSystCorrelations;
  };
}
