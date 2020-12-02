#pragma once

#include "CAFAnaCore/CAFAna/Core/Binning.h"

namespace ana
{
  /// Default true-energy bin edges
  Binning TrueEnergyBins();
  /// Default true-energy bin edges
  const Binning kTrueEnergyBins = TrueEnergyBins();

  // 2km/GeV is the empirical largest L/E at Icarus
  const Binning kTrueLOverEBins = Binning::Simple(100, 0, 2);
}
