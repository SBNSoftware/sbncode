#include "OscLib/IOscCalculator.h"

namespace osc
{
  TMD5* IOscCalculatorAdjustable::GetParamsHashDefault(const std::string& txt) const
  {
    TMD5* ret = new TMD5;
    ret->Update((unsigned char*)txt.c_str(), txt.size());
    const int kNumParams = 8;
    double buf[kNumParams] = {fRho, fL, fDmsq21, fDmsq32,
                              fTh12, fTh13, fTh23, fdCP};
    ret->Update((unsigned char*)buf, sizeof(double)*kNumParams);
    ret->Final();
    return ret;
  }
}
