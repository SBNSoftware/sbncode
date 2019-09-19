#ifndef OSC_OSCCALCULATORPMNS_H
#define OSC_OSCCALCULATORPMNS_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// \file   OscCalculatorPMNS.h                                          //
//                                                                      //
// \brief  Adapt the PMNS calculator to standard interface              //
// \author <bckhouse@caltech.edu>					//
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "OscLib/IOscCalculator.h"
#include "OscLib/PMNS.h"

namespace osc
{
  /// Adapt the \ref PMNS calculator to standard interface
  class OscCalculatorPMNS: public IOscCalculatorAdjustable
  {
  public:
    OscCalculatorPMNS();
    virtual ~OscCalculatorPMNS();

    virtual IOscCalculatorAdjustable* Copy() const override;

    virtual double P(int flavBefore, int flavAfter, double E) override;

    virtual void SetL     (double L     ) override {fPropDirty = true; fL      = L;}
    virtual void SetRho   (double rho   ) override {fPropDirty = true; fRho    = rho;}
    virtual void SetDmsq21(double dmsq21) override {fDmDirty   = true; fDmsq21 = dmsq21;}
    virtual void SetDmsq32(double dmsq32) override {fDmDirty   = true; fDmsq32 = dmsq32;}
    virtual void SetTh12  (double th12  ) override {fMixDirty  = true; fTh12   = th12;}
    virtual void SetTh13  (double th13  ) override {fMixDirty  = true; fTh13   = th13;}
    virtual void SetTh23  (double th23  ) override {fMixDirty  = true; fTh23   = th23;}
    virtual void SetdCP   (double dCP   ) override {fMixDirty  = true; fdCP    = dCP;}

    virtual TMD5* GetParamsHash() const override
    {
      return IOscCalculatorAdjustable::GetParamsHashDefault("PMNS");
    }
  protected:
    PMNS fPMNS;

    bool fMixDirty;
    bool fDmDirty;
    bool fPropDirty;
    double fPrevE;
    int fPrevAnti;
  };

} // namespace

#endif
