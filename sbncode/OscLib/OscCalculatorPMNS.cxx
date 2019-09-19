/////////////////////////////////////////////////////////////////////////////
// \file    OscCalculatorPMNS.cxx 
// \brief   Source file for PMNS oscillation calculation
// \version $Id: OscCalculatorPMNS.cxx,v 1.4 2012-09-25 04:51:35 gsdavies Exp $
// \author  
/////////////////////////////////////////////////////////////////////////////
#include "OscLib/OscCalculatorPMNS.h"

#include <cassert>
#include <cstdlib>

namespace osc
{
  OscCalculatorPMNS::OscCalculatorPMNS()
    : fMixDirty(true), fDmDirty(true), fPropDirty(true), fPrevAnti(0)
  {
  }

  OscCalculatorPMNS::~OscCalculatorPMNS()
  {
  }

  IOscCalculatorAdjustable* OscCalculatorPMNS::Copy() const
  {
    return new OscCalculatorPMNS(*this);
  }

  double OscCalculatorPMNS::P(int flavBefore, int flavAfter, double E)
  {
    const int anti = (flavBefore > 0) ? +1 : -1;
    assert(flavAfter/anti > 0);
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
} // namespace
