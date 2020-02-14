#include "CAFAna/Experiment/GaussianConstraint.h"

#include "CAFAna/Vars/FitVars.h"

#include "CAFAna/Core/MathUtil.h"

namespace ana
{
  //----------------------------------------------------------------------
  double GaussianConstraint::ChiSq(osc::IOscCalculatorAdjustable* osc,
				   const SystShifts& /*syst*/) const
  {
    return util::sqr((fVar->GetValue(osc)-fMean)/fSigma);
  }

}
