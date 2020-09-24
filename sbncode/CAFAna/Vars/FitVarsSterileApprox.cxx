#include "CAFAna/Vars/FitVarsSterileApprox.h"

#include "CAFAna/Core/OscCalcSterileApprox.h"

namespace ana
{
  // --------------------------------------------------------------------------
  double FitDmSqSterile::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return DowncastToSterileApprox(osc)->GetDmsq();
  }

  // --------------------------------------------------------------------------
  void FitDmSqSterile::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    DowncastToSterileApprox(osc)->SetDmsq(Clamp(val));
  }

  // --------------------------------------------------------------------------
  double FitSinSq2ThetaMuMu::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return DowncastToSterileApprox(osc)->GetSinSq2ThetaMuMu();
  }

  // --------------------------------------------------------------------------
  void FitSinSq2ThetaMuMu::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    DowncastToSterileApprox(osc)->SetSinSq2ThetaMuMu(Clamp(val));
  }

  // --------------------------------------------------------------------------
  double FitSinSq2ThetaMuE::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return DowncastToSterileApprox(osc)->GetSinSq2ThetaMuE();
  }

  // --------------------------------------------------------------------------
  void FitSinSq2ThetaMuE::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    DowncastToSterileApprox(osc)->SetSinSq2ThetaMuE(Clamp(val));
  }
}
