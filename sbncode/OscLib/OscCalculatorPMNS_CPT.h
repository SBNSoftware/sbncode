#ifndef OSC_OSCCALCULATORPMNSCPT_H
#define OSC_OSCCALCULATORPMNSCPT_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// \file   OscCalculatorPMNS_CPT.h                                           //
//                                                                      //
// \brief  Adapt the PMNS calculator to standard interface and include  //
//         neutrino and anti neutrino oscillations seperately           //
// \author <m.tamsett@sussex.ac.uk>                                     //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "OscLib/IOscCalculator.h"
#include "OscLib/PMNS.h"

namespace osc
{
  /// Adapt the \ref PMNS calculator to standard interface and include neutrino and anti neutrino oscillations seperately
  class OscCalculatorPMNS_CPT: public IOscCalculatorAdjustable
  {
  public:
    OscCalculatorPMNS_CPT();
    virtual ~OscCalculatorPMNS_CPT();

    virtual IOscCalculatorAdjustable* Copy() const override;

    virtual double P(int flavBefore, int flavAfter, double E) override;

    virtual void SetL     (double L     ) override {fPropDirty = true; fL      = L;}
    virtual void SetRho   (double rho   ) override {fPropDirty = true; fRho    = rho;}
    // Neutrino parameters
    virtual void SetDmsq21(double dmsq21) override {fDmDirty   = true; fDmsq21 = dmsq21;}
    virtual void SetDmsq32(double dmsq32) override {fDmDirty   = true; fDmsq32 = dmsq32;}
    virtual void SetTh12  (double th12  ) override {fMixDirty  = true; fTh12   = th12;}
    virtual void SetTh13  (double th13  ) override {fMixDirty  = true; fTh13   = th13;}
    virtual void SetTh23  (double th23  ) override {fMixDirty  = true; fTh23   = th23;}
    virtual void SetdCP   (double dCP   ) override {fMixDirty  = true; fdCP    = dCP;}
    // Anti neutrino parameters
    virtual void SetDmsq21Bar(double dmsq21_bar){fDmDirty_bar   = true; fDmsq21_bar = dmsq21_bar;}
    virtual void SetDmsq32Bar(double dmsq32_bar){fDmDirty_bar   = true; fDmsq32_bar = dmsq32_bar;}
    virtual void SetTh12Bar  (double th12_bar  ){fMixDirty_bar  = true; fTh12_bar   = th12_bar;}
    virtual void SetTh13Bar  (double th13_bar  ){fMixDirty_bar  = true; fTh13_bar   = th13_bar;}
    virtual void SetTh23Bar  (double th23_bar  ){fMixDirty_bar  = true; fTh23_bar   = th23_bar;}
    virtual void SetdCPBar   (double dCP_bar   ){fMixDirty_bar  = true; fdCP_bar    = dCP_bar;}
    virtual double GetDmsq21Bar() const { return fDmsq21_bar ; }
    virtual double GetDmsq32Bar() const { return fDmsq32_bar ; }
    virtual double GetTh12Bar  () const { return fTh12_bar   ; }
    virtual double GetTh13Bar  () const { return fTh13_bar   ; }
    virtual double GetTh23Bar  () const { return fTh23_bar   ; }
    virtual double GetdCPBar   () const { return fdCP_bar    ; }

    virtual TMD5* GetParamsHash() const override
    {
      return IOscCalculatorAdjustable::GetParamsHashDefault("PMNS");
    }
    TMD5* GetParamsHashDefaultBar() const;
  protected:
    PMNS fPMNS;
    bool fMixDirty;
    bool fDmDirty;
    bool fPropDirty;
    double fPrevE;
    int fPrevAnti;
    
    double fDmsq21_bar;
    double fDmsq32_bar;
    double fTh12_bar;
    double fTh13_bar;
    double fTh23_bar;
    double fdCP_bar;
    PMNS fPMNS_bar;
    bool fMixDirty_bar;
    bool fDmDirty_bar;
    bool fPropDirty_bar;
    double fPrevE_bar;
    int fPrevAnti_bar;
  };

} // namespace

#endif
