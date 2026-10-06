#ifndef SBNCODE_BEAMSPILLRETRIEVER_GETFOM_H
#define SBNCODE_BEAMSPILLRETRIEVER_GETFOM_H

/**
 * @file   sbncode/BeamSpillInfoRetriever/getFOM.h
 * @brief  Beam quality figures of merit.
 * @author Max Dubnowski (maxdub@upenn.sas.edu)
 */

#include "sbnobj/Common/POTAccounting/BNBSpillInfo.h"

#include <tuple>



namespace sbn
{
  /// Bits of the beam-quality status word (same values as the analysis-level
  /// `fom_status` of bnb_fom.py / fom_fixups.py).
  namespace fomstatus {
    constexpr unsigned int NoTOR              = 1;      ///< no valid, positive TOR860 or TOR875 (no beam)
    constexpr unsigned int NoHBPM             = 2;      ///< horizontal position/angle cannot be formed
    constexpr unsigned int NoVBPM             = 4;      ///< vertical position/angle cannot be formed
    constexpr unsigned int NoMWWidth          = 64;     ///< no usable multiwire width: nominal width
    constexpr unsigned int NeighborFilled     = 128;    ///< position/angle from the adjacent spills
    constexpr unsigned int TOR875Used         = 256;    ///< TOR860 drop-out: intensity from TOR875
    constexpr unsigned int MWEmpty            = 512;    ///< database multiwire width is an empty chamber
    constexpr unsigned int WidthFromNeighbors = 1024;   ///< width from the adjacent spills
    constexpr unsigned int FallbackRemoved    = 2048;   ///< only the non-primary BPM read: not used
    constexpr unsigned int BurstFill          = 16384;  ///< 875-station drop-out filled from both sides
    constexpr unsigned int HPTG2Primary       = 32768;  ///< horizontal projection through HPTG2
  }

  /// TOR860 below this fraction of TOR875 is a TOR860 drop-out.
  constexpr double TOR860DropoutFraction = 0.5;

  /// Beam at the centre of the target, as used by the figure of merit.
  struct BNBBeamState {
    enum WidthSource_t: int { FitM876 = 1, FitM875 = 2, DatabaseM876 = 3, DatabaseM875 = 4, Nominal = 5 };
    double tor  = -1.;    ///< intensity used [protons]
    double hpos = -999.;  ///< horizontal position [mm]
    double hang = -999.;  ///< horizontal angle [mrad]
    double vpos = -999.;  ///< vertical position [mm]
    double vang = -999.;  ///< vertical angle [mrad]
    double sx   = -999.;  ///< horizontal width at the target [mm] (-999: nominal)
    double sy   = -999.;  ///< vertical width at the target [mm] (-999: nominal)
    int widthSource = Nominal;
    unsigned int status = 0; ///< bits from `sbn::fomstatus`
    bool hasFOM() const
      { return (status & (fomstatus::NoTOR | fomstatus::NoHBPM | fomstatus::NoVBPM)) == 0; }
  };

  /// Whether HPTG2 (instead of HPTG1) is the primary horizontal target BPM at
  /// this time: the one autotune was calibrated on most recently.
  bool hptg2IsPrimary(unsigned long spill_time_s);

  /// Intensity for the FOM [protons]: TOR860, or TOR875 if TOR860 is missing or
  /// dropped out (`*tor875used` is set in the latter case); -1 if none.
  double beamIntensity(BNBSpillInfo const& spill, bool* tor875used = nullptr);

  /// Position, angle, width and status of the beam at the target for one spill.
  BNBBeamState getBNBBeamState(BNBSpillInfo const& spill);

  /// FOM of a beam state (-999 if it has none); `nominalWidth` ignores its width.
  double computeFOM(BNBBeamState const& state, bool nominalWidth = false);

  /**
   * @brief Returns a Figure of Merit on BNB beam quality.
   *
   * The figure of merit is described in [SBN DocDB 41901](https://sbn-docdb.fnal.gov/cgi-bin/sso/ShowDocument?docid=41901).
   * Inputs the BNBSpillInfo and returns the BNB Quality Metric called FOM, derived from MicroBooNE's FOM
   *
   * @return { FOM with fitted multiwire width, FOM with database ("pre-fit")
   *           width, FOM with nominal width }
   *
   * Each value is in [0, 1], or one of these codes (same for all three):
   *  * `-1`: no valid intensity (TOR860 and TOR875 missing or <= 0, i.e. no beam)
   *  * `2`:  horizontal position/angle cannot be formed (missing BPM or BPM offset)
   *  * `3`:  vertical position/angle cannot be formed (missing BPM or BPM offset)
   * The first two values are additionally `-999` when no usable width is found.
   * A device value of `-999` is treated as missing.
   *
   * Horizontal: HP875 and the primary target BPM (`hptg2IsPrimary()`);
   * vertical: VP875 and VP873. The other target BPMs are not used as fallbacks.
   * Widths: M876 then M875 profile fits (not the target multiwire), then the
   * database widths unless the chamber was empty (sigma > 4 mm).
   * Intensity: `beamIntensity()`.
   */
  std::tuple<float, float, float> getBNBqualityFOM(BNBSpillInfo const& spill);

  /**
    * @brief Inside the getFOM script, takes the positions and angles of the beam and calculates the BNB FOM
    */
  double calcFOM(double horpos,double horang,double verpos,double verang,double tor,double tgtsx=-999,double tgtsy=-999);
    
  /**  
    * @brief Takes in the centroid and sigma of the beam, along with transfer matrices, and will determine the beam's 
    * 2D gaussian position depending where on the target is being measured. The code "swims" up the target to calculate these
    */
  void swimBNB(const double centroid1[6], const double sigma1[6][6], 
               const double xferc[6][6], const double xfers[6][6],
               double &cx, double& cy, double &sx, double &sy, double &rho);
 
  /**
    * @brief Integrates the 2D modelled gaussian beam overlapping with the target, and returns the fraction outside the beam 
    */
  double func_intbivar(const double cx, const double cy, const double sx, const double sy, const double rho );
 
  /**
    * @brief Inputs the MWR Data and determines the centroid, sigma, and chi2 value of a gaussian fit of the beam
    * @return whether the fit succeeded (valid status, positive degrees of freedom)
    */
  bool processBNBprofile(const double* mwdata, double &x, double& sx, double& chi2);
  
}
#endif
