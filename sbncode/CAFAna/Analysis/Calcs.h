#pragma once

namespace osc
{
  template<class T> class _IOscCalcAdjustable;
  typedef _IOscCalcAdjustable<double> IOscCalcAdjustable;
  class OscCalcSterile;
}

namespace ana
{
  /// Reset calc to default assumptions for all parameters
  void ResetOscCalcToDefault(osc::IOscCalcAdjustable* calc);
  void ResetOscCalcToDefaultIH(osc::IOscCalcAdjustable* calc);

  /// Create a new calc with default assumptions for all parameters
  osc::IOscCalcAdjustable* DefaultOscCalc();
  osc::IOscCalcAdjustable* DefaultOscCalcIH();

  /// Reset calc to default assumptions for all parameters
  void ResetSterileCalcToDefault(osc::OscCalcSterile* calc);

  /// Create a sterile calc with default assumptions for all parameters
  osc::OscCalcSterile* DefaultSterileCalc(int nflavors);
}
