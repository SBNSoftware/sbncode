#ifndef SBNCODE_BEAMSPILLRETRIEVER_BNBFOMFILL_H
#define SBNCODE_BEAMSPILLRETRIEVER_BNBFOMFILL_H

/**
 * @file   sbncode/BeamSpillInfoRetriever/BNBFOMFill.h
 * @brief  BNB figure of merit for spills with missing inputs, from their neighbours.
 *
 * These use the other spills of the same subrun, so they run at the end of
 * the subrun in the retriever modules. Validated with closure tests on
 * ICARUS Run 2 and SBND Run 1 (hide what the fill replaces in measured spills,
 * rebuild the FOM, count bad spills, true FOM <= 0.98, that get accepted).
 */

#include "sbnobj/Common/POTAccounting/BNBSpillInfo.h"

#include <vector>

namespace sbn {

  struct BNBFOMFillConfig {

    double minTor = 1e11; ///< spills below this intensity [protons] are not filled

    // ---- empty multiwire / neighbour width
    /// Spills whose database width is an empty chamber take the width of the
    /// adjacent spills when those have a measured width agreeing within
    /// `nbMaxDSig` (otherwise they keep the nominal width).
    bool neighbourWidth = true;

    // ---- neighbour fill: BPM reading missing, adjacent spills stable
    bool neighbourFill = true;
    double nbMaxGap  = 0.1;   ///< each adjacent spill at most this far [s] (next 15 Hz spill)
    double nbMaxDTor = 0.03;  ///< spill intensity within this fraction of the neighbours' mean
    double nbMaxDPos = 1.0;   ///< neighbours' positions agree [mm]
    double nbMaxDAng = 0.4;   ///< neighbours' angles agree [mrad]
    double nbMaxDSig = 0.1;   ///< neighbours' widths agree [mm]
    double nbMaxDFOM = 0.001; ///< neighbours' FOMs agree
    double nbMinFOM  = -1.;   ///< both neighbours' FOM above this (< 0: no requirement; ICARUS: 0.998)

    // ---- 875-station drop-outs (HP875 and VP875 missing, target BPMs read)
    bool burstFill = false;     ///< ICARUS
    double burstMaxDt   = 60.;    ///< measured spill on each side within this [s]
    double burstMinFOM  = 0.999;  ///< both of them above this FOM
    double burstMaxDAng = 0.2;    ///< their angles agree [mrad]
    double burstMaxDTgt = 0.5;    ///< own HPTG1/VPTG2 within this of their mean [mm]
    double burstAccept  = 0.995;  ///< filled only if the estimated FOM is above this

  }; // BNBFOMFillConfig

  /**
   * @brief Improves the FOM of the spills of one subrun using their neighbours.
   * @param spills all the spills of the subrun (any order); FOMs are updated
   * @param config which fills to apply and their thresholds
   * @return the status word (`sbn::fomstatus` bits) of each spill, same order
   *
   * Spills that are not filled keep the result of `getBNBqualityFOM()`.
   * A filled spill gets its FOM in `FOM` (measured width) or in
   * `NoMultiWireFOM` (nominal width), the others set to -999; its BPM
   * readings are not changed, so a filled spill still shows missing readings.
   *
   * 1. neighbour width (bit `WidthFromNeighbors`)
   * 2. neighbour fill (bit `NeighborFilled`): the spill takes the mean position
   *    and angle of the spills just before and after (and their width if it
   *    has none), with its own intensity
   * 3. 875 drop-out fill (bit `BurstFill`): own HPTG1/VPTG2 plus the offset to
   *    the target and the angle of the nearest measured spill on each side;
   *    BPM offsets cancel in the differences
   */
  std::vector<unsigned int> improveBNBqualityFOMs
    (std::vector<BNBSpillInfo>& spills, BNBFOMFillConfig const& config);

} // namespace sbn

#endif // SBNCODE_BEAMSPILLRETRIEVER_BNBFOMFILL_H
