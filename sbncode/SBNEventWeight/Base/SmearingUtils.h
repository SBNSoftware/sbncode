#ifndef _SBN_SMEARINGUTILS_H_
#define _SBN_SMEARINGUTILS_H_

/**
 * Multivariate Gaussian sampling utilities.
 *
 * Original author: J. Zennamo.
 *
 * Note: In larsim EventWeight, these functions are static members of the
 * WeightCalc class.
 */

#include "CLHEP/Random/RandGaussQ.h"
#include "TMatrixD.h"
#include <vector>

namespace sbn {
  namespace evwgh {

/**
 * Apply Gaussian smearing to a set of data.
 *
 * If centralValues is of dimension N, inputCovarianceMatrix needs to be NxN,
 * and each of the returned data sets will be also of dimension N.
 *
 * @param centralValues the values to be smeared
 * @param inputCovarianceMatrix covariance matrix for smearing
 * @param n_multisims number of sets of smeared values to be produced
 * @return a set of n_multisims value sets smeared from the central value
 */
std::vector<std::vector<double> > MultiGaussianSmearing(
  std::vector<double> const& centralValues,
  std::vector<std::vector<double> > const& inputCovarianceMatrix,
  int n_multisims,
  CLHEP::RandGaussQ& GaussRandom);

std::vector<double>  MultiGaussianSmearing(
  std::vector<double> const& centralValue,
  TMatrixD* const& inputCovarianceMatrix,
  std::vector<double> rand);

std::vector<double> MultiGaussianSmearing(
  std::vector<double> const& centralValue,
  TMatrixD* const& LowerTriangleCovarianceMatrix,
  bool isDecomposed,
  std::vector<double> rand);

  }  // namespace evwgh
}  // namespace sbn

#endif  // _SBN_SMEARINGUTILS_H_

