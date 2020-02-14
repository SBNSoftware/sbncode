#pragma once

namespace osc{class IOscCalculatorAdjustable;
              class OscCalculatorSterile;}

namespace ana
{
  /// Reset calculator to default assumptions for all parameters
  void ResetOscCalcToDefault(osc::IOscCalculatorAdjustable* calc);
  void ResetOscCalcToDefaultIH(osc::IOscCalculatorAdjustable* calc);

  /// Create a new calculator with default assumptions for all parameters
  osc::IOscCalculatorAdjustable* DefaultOscCalc();
  osc::IOscCalculatorAdjustable* DefaultOscCalcIH();

  /// Reset calculator to default assumptions for all parameters
  void ResetSterileCalcToDefault(osc::OscCalculatorSterile* calc);

  /// Create a sterile calculator with default assumptions for all parameters
  osc::OscCalculatorSterile* DefaultSterileCalc(int nflavors);
}
