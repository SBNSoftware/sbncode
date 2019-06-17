/////////////////////////////////////////////////////////////////////////////
// \file    OscCalculatorPMNS_CPT.cxx 
// \brief   Source file for PMNS oscillation calculation
// \version $Id: OscCalculatorPMNS_CPT.cxx,v 1.0 2013-11-21 11:17:35 tamsett Exp $
// \author  Matthew Tamsett
/////////////////////////////////////////////////////////////////////////////
#include "OscLib/OscCalculatorPMNS_CPT.h"

#include <cassert>
#include <cstdlib>
#include <iostream>

namespace osc
{
  OscCalculatorPMNS_CPT::OscCalculatorPMNS_CPT()
    : fMixDirty(true), fDmDirty(true), fPropDirty(true), fPrevAnti(0),
      fMixDirty_bar(true), fDmDirty_bar(true), fPropDirty_bar(true), fPrevAnti_bar(0)
  {
  }

  OscCalculatorPMNS_CPT::~OscCalculatorPMNS_CPT()
  {
  }

  IOscCalculatorAdjustable* OscCalculatorPMNS_CPT::Copy() const
  {
    return new OscCalculatorPMNS_CPT(*this);
  }

  double OscCalculatorPMNS_CPT::P(int flavBefore, int flavAfter, double E)
  {
    const int anti = (flavBefore > 0) ? +1 : -1;
    assert(flavAfter/anti > 0);
    // We will now forward to the appropriate calculator depending on anti
    if (anti == -1){        
        if(anti != fPrevAnti_bar) fPropDirty_bar = true;    
        int i = -1, j = -1;
        if(abs(flavBefore) == 12) i = 0;
        if(abs(flavBefore) == 14) i = 1;
        if(abs(flavBefore) == 16) i = 2;
        if(abs(flavAfter) == 12) j = 0;
        if(abs(flavAfter) == 14) j = 1;
        if(abs(flavAfter) == 16) j = 2;
        assert(i >= 0 && j >= 0);
    
        if(fMixDirty_bar){
          fPMNS_bar.SetMix(fTh12_bar, fTh23_bar, fTh13_bar, fdCP_bar);
          fMixDirty_bar = false;
        }
        if(fDmDirty_bar){
          fPMNS_bar.SetDeltaMsqrs(fDmsq21_bar, fDmsq32_bar);
          fDmDirty_bar = false;
        }
    
        if(fPropDirty_bar || E != fPrevE_bar){
          fPMNS_bar.Reset();
          // Assume Z/A=0.5
          const double Ne = fRho/2;
          fPMNS_bar.PropMatter(fL, E, Ne, anti);
    
          fPropDirty_bar = false;
          fPrevE_bar = E;
          fPrevAnti_bar = anti;
        }
    
        return fPMNS_bar.P(i, j);
    } else {
        if(anti != fPrevAnti) fPropDirty = true;

        int i = -1, j = -1;
        if(abs(flavBefore) == 12) i = 0;
        if(abs(flavBefore) == 14) i = 1;
        if(abs(flavBefore) == 16) i = 2;
        if(abs(flavAfter) == 12) j = 0;
        if(abs(flavAfter) == 14) j = 1;
        if(abs(flavAfter) == 16) j = 2;
        assert(i >= 0 && j >= 0);
    
        if(fMixDirty){
          fPMNS.SetMix(fTh12, fTh23, fTh13, fdCP);
          fMixDirty = false;
        }
        if(fDmDirty){
          fPMNS.SetDeltaMsqrs(fDmsq21, fDmsq32);
          fDmDirty = false;
        }
    
        if(fPropDirty || E != fPrevE){
          fPMNS.Reset();
          // Assume Z/A=0.5
          const double Ne = fRho/2;
          fPMNS.PropMatter(fL, E, Ne, anti);
    
          fPropDirty = false;
          fPrevE = E;
          fPrevAnti = anti;
        }
    
        return fPMNS.P(i, j);
    }
  }
  TMD5* OscCalculatorPMNS_CPT::GetParamsHashDefaultBar() const
  {
    const std::string& txt = "PMNSBar";
    TMD5* ret = new TMD5;
    ret->Update((unsigned char*)txt.c_str(), txt.size());
    const int kNumParams = 8;
    double buf[kNumParams] = {fRho, fL, fDmsq21_bar, fDmsq32_bar,
                              fTh12_bar, fTh13_bar, fTh23_bar, fdCP_bar};
    ret->Update((unsigned char*)buf, sizeof(double)*kNumParams);
    ret->Final();
    return ret;
  }

} // namespace
