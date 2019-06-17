#include "CAFAna/Vars/FitVars.h"

#include "OscLib/func/IOscCalculator.h"
#include "Utilities/func/MathUtil.h"

#include <cassert>
#include <cmath>

namespace ana
{
  //----------------------------------------------------------------------
  double FitTheta13::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return osc->GetTh13();
  }

  //----------------------------------------------------------------------
  void FitTheta13::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc->SetTh13(val);
  }

  //----------------------------------------------------------------------
  double FitSinSq2Theta13::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return util::sqr(sin(2*osc->GetTh13()));
  }

  //----------------------------------------------------------------------
  void FitSinSq2Theta13::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc->SetTh13(asin(sqrt(Clamp(val)))/2);
  }

  //----------------------------------------------------------------------
  double FitDeltaInPiUnits::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    double ret = osc->GetdCP()/M_PI;

    // convert to the range 0-2
    long long int a = ret/2+1;
    ret -= 2*a;
    // Instead of figuring out all the rounding just nudge the last little bit
    while(ret < 0) ret += 2;
    while(ret > 2) ret -= 2;

    return ret;
  }

  //----------------------------------------------------------------------
  void FitDeltaInPiUnits::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc->SetdCP(M_PI*val);
  }
  //----------------------------------------------------------------------
  double FitTheta23::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return osc->GetTh23();
  }

  //----------------------------------------------------------------------
  void FitTheta23::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc->SetTh23(val);
  }

  //----------------------------------------------------------------------

  //----------------------------------------------------------------------
  double FitSinSqTheta23::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return util::sqr(sin(osc->GetTh23()));
  }

  //----------------------------------------------------------------------
  void FitSinSqTheta23::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc->SetTh23(asin(sqrt(Clamp(val))));
  }

  //----------------------------------------------------------------------
  double FitSinSq2Theta23::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return util::sqr(sin(2*osc->GetTh23()));
  }

  //----------------------------------------------------------------------
  void FitSinSq2Theta23::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc->SetTh23(asin(sqrt(Clamp(val)))/2);
  }

  //----------------------------------------------------------------------
  double FitDmSq32::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return osc->GetDmsq32();
  }

  //----------------------------------------------------------------------
  void FitDmSq32::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc->SetDmsq32(Clamp(val));
  }

  //----------------------------------------------------------------------
  double FitDmSq32Scaled::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return osc->GetDmsq32()*1000.0;
  }

  //----------------------------------------------------------------------
  void FitDmSq32Scaled::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc->SetDmsq32(Clamp(val/1000.0));
  }

  //----------------------------------------------------------------------
  double FitTanSqTheta12::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return util::sqr(tan(osc->GetTh12()));
  }

  //----------------------------------------------------------------------
  void FitTanSqTheta12::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc->SetTh12(atan(sqrt(Clamp(val))));
  }

  //----------------------------------------------------------------------
  double FitSinSq2Theta12::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return util::sqr(sin(2*osc->GetTh12()));
  }

  //----------------------------------------------------------------------
  void FitSinSq2Theta12::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc->SetTh12(asin(sqrt(Clamp(val)))/2);
  }

  //----------------------------------------------------------------------
  double FitDmSq21::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return osc->GetDmsq21();
  }

  //----------------------------------------------------------------------
  void FitDmSq21::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc->SetDmsq21(Clamp(val));
  }
  //----------------------------------------------------------------------
  double FitRho::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return osc->GetRho();
  }

  //----------------------------------------------------------------------
  void FitRho::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc->SetRho(Clamp(val));
  }

} // namespace
