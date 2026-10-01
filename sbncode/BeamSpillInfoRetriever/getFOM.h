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
