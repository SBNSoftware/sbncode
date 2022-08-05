/**
 * \file TruncMean.h
 *
 * \ingroup 3DMichel
 * 
 * \brief Class def header for a class TruncMean
 *
 * @author david caratelli [davidc@fnal.gov]
 * Written 08/02/2018.
 *
 */

#ifndef TRUNCMEAN_H
#define TRUNCMEAN_H

#include <iostream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <climits>
#include <limits>

/**
   \class TruncMean
   The truncated mean class allows to compute the following quantities
   1) the truncated mean profile of an ordered vector of values, such as 
   the charge profile along a particle's track.
   To create such a profile use the function CalcTruncMeanProfile()
   2) Get the truncated mean value of a distribution. This function
   iteratively hones in on the truncated mean of a distribution by
   updating the mean and cutting the tails after each iteration.
   For this functionality use CalcIterativeTruncMean()
   doxygen documentation!
*/

namespace single_photon
{
  static const double kINVALID_FLOAT = std::numeric_limits<double>::max();

  class TruncMean{

    public:

      /// Default constructor
      TruncMean(){}

      /// Default destructor
      ~TruncMean(){}

      /**
        @brief Given residual range and dq vectors return truncated local dq.
        Input vectors are assumed to be match pair-wise (nth entry in rr_v
        corresponds to nth entry in dq_v vector).
        Input rr_v values are also assumed to be ordered: monotonically increasing
        or decreasing.
        For every dq value a truncated linear dq value is calculated as follows:
        0) all dq values within a rr range set by the class variable _rad are selected.
        1) the median and rms of these values is calculated.
        2) the subset of local dq values within the range [median-rms, median+rms] is selected.
        3) the resulting local truncated dq is the average of this truncated subset.
        @input std::vector<double> rr_v -> vector of x-axis coordinates (i.e. position for track profile)
        @input std::vector<double> dq_v -> vector of measured values for which truncated profile is requested
        (i.e. charge profile of a track)
        @input std::vector<double> dq_trunc_v -> passed by reference -> output stored here
        @input double nsigma -> optional parameter, number of sigma to keep around RMS for TM calculation
        */
      void CalcTruncMeanProfile(const std::vector<double>& rr_v, const std::vector<double>& dq_v,
          std::vector<double>& dq_trunc_v, const double& nsigma = 1);

      /**
        @brief Iteratively calculate the truncated mean of a distribution
        @brief: mean is returned if vecter's size is too small, or reach the max iteration, or median has converged
        @input std::vector<double> v -> vector of values for which truncated mean is asked
        @input size_t nmin -> minimum number of iterations to converge on truncated mean
        @input size_t nmax -> maximum number of iterations to converge on truncated mean
        @input size_t lmin -> minimum number of entries in vector before exiting and returning current value
        @input size_t currentiteration -> current iteration
        @input double convergencelimit -> fractional difference between successive iterations
        under which the iteration is completed, provided nmin iterations have occurred.
        @input nsigma -> number of sigma around the median value to keep when the distribution is trimmed.
        */
      double CalcIterativeTruncMean(std::vector<double> v, const size_t& nmin,
          const size_t& nmax, const size_t& currentiteration,
          const size_t& lmin,
          const double& convergencelimit,
          const double& nsigma, const double& oldmed = kINVALID_FLOAT);

      /**
        @brief Set the smearing radius over which to take hits for truncated mean computaton.
        */
      void setRadius(const double& rad) { _rad = rad; }

    private:

      double Mean  (const std::vector<double>& v);
      double Median(const std::vector<double>& v);
      double RMS   (const std::vector<double>& v);

      /**
        Smearing radius over which charge from neighboring hits is scanned to calculate local
        truncated mean
        */
      double _rad;

  };


}
#endif
/** @} */ // end of doxygen group 
