#include "CAFAna/Core/Binning.h"

namespace ana
{
  //----------------------------------------------------------------------
  Binning TrueEnergyBins()
  {
    // Osc P is ~sin^2(1/E). Difference in prob across a bin is ~dP/dE which
    // goes as 1/E^2 times a trigonometric term depending on the parameters but
    // bounded between -1 and +1. So bins should have width ~1/E^2. E_i~1/i
    // achieves that.

    // Binning roughly re-optimized for SBN sterile analyses

    const int kNumTrueEnergyBins = 80;

    // N+1 bin low edges
    std::vector<double> edges(kNumTrueEnergyBins+1);

    const double Emin = .2; // 200 MeV

    // How many edges to generate. Allow room for 0-Emin bin
    const double N = kNumTrueEnergyBins-1;
    const double A = N*Emin;

    edges[0] = 0;

    for(int i = 1; i <= N; ++i){
      edges[kNumTrueEnergyBins-i] = A/i;
    }

    edges[kNumTrueEnergyBins] = 20; // Replace the infinity that would be here

    return Binning::Custom(edges);
  }


  //----------------------------------------------------------------------
  Binning TrueLOverTrueEBins()
  {
    // constant binnig

    const int kNumTrueLOverTrueEBins = 2000;
    const double klow = 0.0;
    const double khigh = 5.0;

    return Binning::Simple(kNumTrueLOverTrueEBins, klow, khigh);
  }
}
