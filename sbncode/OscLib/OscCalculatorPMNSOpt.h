#ifndef OSC_OSCCALCULATORPMNSOPT_H
#define OSC_OSCCALCULATORPMNSOPT_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// \file OscCalculatorPMNSOpt.h                                         //
//                                                                      //
// Adapt the PMNSOpt calculator to standard interface                   //
// <bckhouse@caltech.edu>						//
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "IOscCalculator.h"
#include "PMNSOpt.h"

#include <map>

namespace osc
{
  /// \brief Optimized version of \ref OscCalculatorPMNS
  ///
  /// Adapt the \ref PMNSOpt calculator to standard interface
  class OscCalculatorPMNSOpt: public IOscCalculatorAdjustable
  {
  public:
    OscCalculatorPMNSOpt();
    virtual ~OscCalculatorPMNSOpt();

    virtual IOscCalculatorAdjustable* Copy() const override;

    virtual double P(int flavBefore, int flavAfter, double E) override;

    virtual void SetL     (double L     ) override {++fLRIdx;  fL      = L;}
    virtual void SetRho   (double rho   ) override {++fLRIdx;  fRho    = rho;}
    virtual void SetDmsq21(double dmsq21) override {++fDmIdx;  fDmsq21 = dmsq21;}
    virtual void SetDmsq32(double dmsq32) override {++fDmIdx;  fDmsq32 = dmsq32;}
    virtual void SetTh12  (double th12  ) override {++fMixIdx; fTh12   = th12;}
    virtual void SetTh13  (double th13  ) override {++fMixIdx; fTh13   = th13;}
    virtual void SetTh23  (double th23  ) override {++fMixIdx; fTh23   = th23;}
    virtual void SetdCP   (double dCP   ) override {++fMixIdx; fdCP    = dCP;}

    virtual TMD5* GetParamsHash() const override
    {
      return IOscCalculatorAdjustable::GetParamsHashDefault("PMNSOpt");
    }
  protected:
    // How many times the mixing parameters and splittings have been set
    long fMixIdx;
    long fDmIdx;
    long fLRIdx;

    struct Val_t
    {
      Val_t() : mixIdx(-1), dmIdx(-1), lrIdx(-1), pmns(0) {}

      // How many times the mixing parameters and splittings had been set when
      // 'pmns' was last updated. If too small then 'pmns' must be updated
      // before use.
      long mixIdx;
      long dmIdx;
      long lrIdx;
      double P[3][3]; ///< Cache of oscillation probabilities
      PMNSOpt* pmns;  ///< The calculator itself
    };

    std::map<double, Val_t> fPMNSOpt[2]; // [anti][E]
  };

} // namespace

#endif
