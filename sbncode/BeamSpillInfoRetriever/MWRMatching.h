#ifndef SBNCODE_BEAMSPILLINFORETRIEVER_MWRMATCHING_H
#define SBNCODE_BEAMSPILLINFORETRIEVER_MWRMATCHING_H

/**
 * @file   sbncode/BeamSpillInfoRetriever/MWRMatching.h
 * @brief  Association of multiwire chamber readings to BNB spills.
 *
 * Factored out of the ICARUS/SBND BNB retriever modules, which carried three
 * copies of the same loop. Kept free of art dependencies so it can be tested
 * standalone.
 */

#include <cstddef>
#include <vector>

namespace sbn::pot {

  /**
   * @brief Finds, for each multiwire device, the reading belonging to spill `i`.
   * @param MWR_times    reading times per device [s] (already corrected for the
   *                     MWR-to-toroid delay)
   * @param spill_times  toroid spill times [s]
   * @param i            index in `spill_times` of the spill to match
   * @param windowLow    spills with time <= windowLow are not considered as
   *                     competitors
   * @param windowHigh   spills with time > windowHigh are not considered as
   *                     competitors
   * @param maxTimeDiff  largest accepted |t_MWR - t_spill| [s]; non-positive
   *                     disables the requirement
   * @return one index per device into `MWR_times[dev]`, or `-1` if no reading
   *         can be associated with this spill
   *
   * A reading is a candidate for spill `i` only if no other spill in the
   * window is closer to it in time; among candidates the closest one wins.
   *
   * Differences from the original inline code:
   *  - no reading -> `-1` (it used to silently fall back to index 0, attaching
   *    an arbitrary, possibly far away, profile to the spill);
   *  - the best candidate is rejected if it is more than `maxTimeDiff` away
   *    (there was no limit: profiles tens of seconds away were being used).
   */
  std::vector<int> matchMWRToSpill(std::vector<std::vector<double>> const& MWR_times,
                                   std::vector<double> const& spill_times,
                                   std::size_t i,
                                   double windowLow, double windowHigh,
                                   double maxTimeDiff);

} // namespace sbn::pot

#endif
