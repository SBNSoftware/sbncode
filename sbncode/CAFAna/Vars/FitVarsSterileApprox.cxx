#include "CAFAna/Vars/FitVarsSterileApprox.h"

#include "CAFAna/Core/OscCalcSterileApprox.h"

namespace ana
{
  // --------------------------------------------------------------------------
  double FitDmSqSterile::GetValue(const osc::IOscCalcAdjustable* osc) const
  {
    return DowncastToSterileApprox(osc)->GetDmsq();
  }

  // --------------------------------------------------------------------------
  void FitDmSqSterile::SetValue(osc::IOscCalcAdjustable* osc, double val) const
  {
    DowncastToSterileApprox(osc)->SetDmsq(Clamp(val));
  }

  // --------------------------------------------------------------------------
  double FitSinSq2ThetaMuMu::GetValue(const osc::IOscCalcAdjustable* osc) const
  {
    return DowncastToSterileApprox(osc)->GetSinSq2ThetaMuMu();
  }

  // --------------------------------------------------------------------------
  void FitSinSq2ThetaMuMu::SetValue(osc::IOscCalcAdjustable* osc, double val) const
  {
    DowncastToSterileApprox(osc)->SetSinSq2ThetaMuMu(Clamp(val));
  }

  // --------------------------------------------------------------------------
  double FitSinSq2ThetaMuE::GetValue(const osc::IOscCalcAdjustable* osc) const
  {
    return DowncastToSterileApprox(osc)->GetSinSq2ThetaMuE();
  }

  // --------------------------------------------------------------------------
  void FitSinSq2ThetaMuE::SetValue(osc::IOscCalcAdjustable* osc, double val) const
  {
    DowncastToSterileApprox(osc)->SetSinSq2ThetaMuE(Clamp(val));
  }
}
