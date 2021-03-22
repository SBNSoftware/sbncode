#include "CLHEP/Random/RandGaussQ.h"
#include "TMatrixD.h"
#include "TDecompChol.h"
#include "canvas/Utilities/Exception.h"
#include "sbncode/SBNEventWeight/Base/WeightCalc.h"

namespace sbn {
  namespace evwgh {

std::vector<std::vector<double> > MultiGaussianSmearing(
    std::vector<double> const& centralValue,
    std::vector< std::vector<double> > const& inputCovarianceMatrix,
    int n_multisims,
    CLHEP::RandGaussQ& GaussRandom) {
  std::vector<std::vector<double> > setOfSmearedCentralValues;

  // Check that covarianceMatrix is of compatible dimension with the central values
  unsigned int covarianceMatrix_dim = centralValue.size();

  if (inputCovarianceMatrix.size() != covarianceMatrix_dim) {
    throw art::Exception(art::errors::Configuration)
      << "inputCovarianceMatrix rows " << inputCovarianceMatrix.size()
      << " not equal to entries of centralValue[]: " << covarianceMatrix_dim;
  }

  for (auto const & row : inputCovarianceMatrix) {
    if (row.size() != covarianceMatrix_dim) {
      throw art::Exception(art::errors::Configuration)
        << "inputCovarianceMatrix columns " << row.size()
        << " not equal to entries of centralValue[]: " << covarianceMatrix_dim;
     }
  }

  // Change covariance matrix into a TMatrixD object
  int dim = int(inputCovarianceMatrix.size());
  TMatrixD covarianceMatrix(dim, dim);
  for (size_t i=0; i<inputCovarianceMatrix.size(); i++) {
    for (size_t j=0; j<inputCovarianceMatrix[i].size(); j++) {
      covarianceMatrix[i][j] = inputCovarianceMatrix[i][j];
    }
  }

  // Perform Choleskey Decomposition
  TDecompChol dc = TDecompChol(covarianceMatrix);
  if(!dc.Decompose()) {
    throw art::Exception(art::errors::StdException)
      << "Cannot decompose covariance matrix to begin smearing.";
    return setOfSmearedCentralValues;
  }

  // Get upper triangular matrix. This maintains the relations in the
  // covariance matrix, but simplifies the structure.
  TMatrixD U = dc.GetU();

  // for every multisim
  for (int n=0; n<n_multisims; n++) {
    // Get a gaussian random number for every central value element
    int dim = centralValue.size();

    std::vector<double> rands(dim);
    GaussRandom.fireArray(dim, rands.data());

    // Compute the smeared central values
    std::vector<double> smearedCentralValues;
    for(int col=0; col<dim; col++) {
      // Find the weight of each col of the upper triangular cov. matrix
      double weightFromU = 0.;
      for (int row=0; row<col+1; row++) {
        weightFromU += U(row,col) * rands[row];
      }

      // multiply this weight by each col of the central values to obtain
      // the gaussian smeared and constrained central values
      // smearedCentralValues.push_back(weightFromU * centralValue[col]);
      smearedCentralValues.push_back(weightFromU + centralValue[col]);
    }

    // collect the smeared central values into a set
    setOfSmearedCentralValues.push_back(smearedCentralValues);
  }

  return setOfSmearedCentralValues;
}


std::vector<double> MultiGaussianSmearing(
    std::vector<double> const& centralValue,
    TMatrixD* const& inputCovarianceMatrix,
    std::vector< double > rand) {
  std::vector<double> smearedCentralValues;

  // Perform Choleskey Decomposition
  // see http://pdg.lbl.gov/2016/reviews/rpp2016-rev-monte-carlo-techniques.pdf (Page 5)
  TDecompChol dc = TDecompChol(*(inputCovarianceMatrix));
  if (!dc.Decompose()) {
    throw art::Exception(art::errors::StdException)
      << "Cannot decompose covariance matrix to begin smearing.";
    return smearedCentralValues;
  }

  // Get upper triangular matrix. This maintains the relations in the
  // covariance matrix, but simplifies the structure.
  TMatrixD U = dc.GetU();

  for (size_t col=0; col<centralValue.size(); col++) {
    // Find the weight of each col of the upper triangular cov. matrix
    double weightFromU = 0;
    for (unsigned int row=0; row<col+1; row++) {
      weightFromU += U(row,col) * rand[row];
    }

    // multiply this weight by each col of the central values to obtain
    // the gaussian smeared and constrained central values
    // smearedCentralValues.push_back(weightFromU * centralValue[col]);
    smearedCentralValues.push_back(weightFromU + centralValue[col]);
  }

  return smearedCentralValues;
}


std::vector<double> MultiGaussianSmearing(
    std::vector<double> const& centralValue,
    TMatrixD* const& LowerTriangleCovarianceMatrix,
     bool isDecomposed,
     std::vector< double > rand) {
  std::vector<double> smearedCentralValues;

  if (!isDecomposed) {
    throw art::Exception(art::errors::StdException)
      << "Must supply the decomposed lower triangular covariance matrix.";
    return smearedCentralValues;
  }

  for(unsigned int col = 0; col < centralValue.size(); ++col) {
    // Find the weight of each col of the upper triangular cov. matrix
    double weightFromU = 0;
    for(size_t row=0; row<col+1; row++) {
      weightFromU += LowerTriangleCovarianceMatrix[0][row][col] * rand[row];
    }

    // multiply this weight by each col of the central values to obtain
    // the gaussian smeared and constrained central values
    // smearedCentralValues.push_back(weightFromU * centralValue[col]);
    smearedCentralValues.push_back(weightFromU + centralValue[col]);
  }

  return smearedCentralValues;
}

  }  // namespace evwgh
}  // namespace sbn

