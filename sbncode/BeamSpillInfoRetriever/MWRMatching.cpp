/**
 * @file   sbncode/BeamSpillInfoRetriever/MWRMatching.cpp
 * @brief  Association of multiwire chamber readings to BNB spills.
 */
#include "sbncode/BeamSpillInfoRetriever/MWRMatching.h"

#include <cmath>
#include <limits>

std::vector<int> sbn::pot::matchMWRToSpill(std::vector<std::vector<double>> const& MWR_times,
                                           std::vector<double> const& spill_times,
                                           std::size_t i,
                                           double windowLow, double windowHigh,
                                           double maxTimeDiff)
{
  std::vector<int> matched(MWR_times.size(), -1);
  if (i >= spill_times.size()) return matched;
  double const t_spill = spill_times[i];

  for (std::size_t dev = 0; dev < MWR_times.size(); ++dev) {
    double Tdiff = std::numeric_limits<double>::max();
    for (std::size_t mwrt = 0; mwrt < MWR_times[dev].size(); ++mwrt) {
      double const t_mwr = MWR_times[dev][mwrt];
      double const d = std::abs(t_mwr - t_spill);
      if (d >= Tdiff) continue;

      // is another spill in the window a better match for this reading?
      bool best_match = true;
      for (std::size_t j = 0; j < spill_times.size(); ++j) {
        if (j == i) continue;
        if (spill_times[j] > windowHigh) continue;
        if (spill_times[j] <= windowLow) continue;
        if (std::abs(t_mwr - spill_times[j]) < d) { best_match = false; break; }
      }
      if (best_match) {
        matched[dev] = static_cast<int>(mwrt);
        Tdiff = d;
      }
    }
    if (matched[dev] >= 0 && maxTimeDiff > 0. && Tdiff > maxTimeDiff) matched[dev] = -1;
  }
  return matched;
}
