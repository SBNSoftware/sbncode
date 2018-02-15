#include "sbncode/EventWeight/WeightCalc.h"

// art libraries
#include "canvas/Utilities/Exception.h"
//#include "art/Utilities/Exception.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// ROOT libraries
#include "TMatrixD.h"
#include "TDecompChol.h"

// CLHEP libraries
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGaussQ.h"

namespace sbncode {
namespace evwgh { 
  std::vector<std::vector<double> > WeightCalc::MultiGaussianSmearing(std::vector<double> const& centralValue,std::vector< std::vector<double> > const& inputCovarianceMatrix,int n_multisims,CLHEP::RandGaussQ& GaussRandom)
  {

    std::vector<std::vector<double> > setOfSmearedCentralValues;

    //Check that covarianceMatrix is of compatible dimension with the central values
    unsigned int covarianceMatrix_dim = centralValue.size();
 
    if (inputCovarianceMatrix.size() != covarianceMatrix_dim)
      {
	throw art::Exception(art::errors::Configuration)
	  << "inputCovarianceMatrix rows " << inputCovarianceMatrix.size()
	  << " not equal to entries of centralValue[]: " << covarianceMatrix_dim;
      }

    for (auto const & row : inputCovarianceMatrix )
      { // Range-for (C++11).
	if(row.size() != covarianceMatrix_dim)
	  {
	    throw art::Exception(art::errors::Configuration)
	      << "inputCovarianceMatrix columns " << row.size()
	      << " not equal to entries of centralValue[]: " << covarianceMatrix_dim;
	  }
      }

    //change covariance matrix into a TMartixD object
    int dim = int(inputCovarianceMatrix.size());
    TMatrixD covarianceMatrix(dim,dim);
    int i = 0;
    for (auto const & row : inputCovarianceMatrix)
      {
	int j = 0;
	for (auto const & element : row)
	  {
	    covarianceMatrix[i][j] = element;
	    ++j;
	  }
    	++i;
      }

    //perform Choleskey Decomposition
    TDecompChol dc = TDecompChol(covarianceMatrix);
    if(!dc.Decompose())
      {
    	throw art::Exception(art::errors::StdException)
	  << "Cannot decompose covariance matrix to begin smearing.";
	return setOfSmearedCentralValues;
      }

    //Get upper triangular matrix. This maintains the relations in the
    //covariance matrix, but simplifies the structure.
    TMatrixD U = dc.GetU();

    //for every multisim
    for(int n = 0; n < n_multisims; ++n)
      {

    	//get a gaussian random number for every central value element
	int dim = centralValue.size();
      
	std::vector<double> rands(dim);
	GaussRandom.fireArray(dim, rands.data());
      
	//compute the smeared central values
	std::vector<double> smearedCentralValues;
	for(int col = 0; col < dim; ++col)
	  {
	    //find the weight of each col of the upper triangular cov. matrix
	    double weightFromU = 0.;
	    for(int row = 0; row < col+1; ++row)
	      {
        	weightFromU += U(row,col)*rands[row];
	      }
	    //multiply this weight by each col of the central values to obtain
	    //the gaussian smeared and constrained central values
	    // smearedCentralValues.push_back(weightFromU * centralValue[col]);
	    smearedCentralValues.push_back(weightFromU + centralValue[col]);
	  }

	//collect the smeared central values into a set
	setOfSmearedCentralValues.push_back(smearedCentralValues);
      }//loop over multisims
    return setOfSmearedCentralValues;
  } // WeightCalc::MultiGaussianSmearing() - gives a vector of sets of smeared parameters

  /////////////
  /// Over load of the above function that only returns a single varied parameter set
  ////////////
  std::vector<double> WeightCalc::MultiGaussianSmearing(std::vector<double> const& centralValue, TMatrixD* const& inputCovarianceMatrix, std::vector< double > rand)
  {

    //compute the smeared central values
    std::vector<double> smearedCentralValues;

    //
    //perform Choleskey Decomposition
    //
    // Best description I have found 
    //     is in the PDG (Monte Carlo techniques, Algorithms, Gaussian distribution)
    //
    //  http://pdg.lbl.gov/2016/reviews/rpp2016-rev-monte-carlo-techniques.pdf (Page 5)
    //
    //
    TDecompChol dc = TDecompChol(*(inputCovarianceMatrix));
    if(!dc.Decompose())
      {
    	throw art::Exception(art::errors::StdException)
	  << "Cannot decompose covariance matrix to begin smearing.";
	return smearedCentralValues;
      }

    //Get upper triangular matrix. This maintains the relations in the
    //covariance matrix, but simplifies the structure.
    TMatrixD U = dc.GetU();

 
    for(unsigned int col = 0; col < centralValue.size(); ++col)
      {
	//find the weight of each col of the upper triangular cov. matrix
	double weightFromU = 0.;
	
	for(unsigned int row = 0; row < col+1; ++row)
	  {
	    weightFromU += U(row,col)*rand[row];
 	  }

	//multiply this weight by each col of the central values to obtain
	//the gaussian smeared and constrained central values
	// smearedCentralValues.push_back(weightFromU * centralValue[col]);
	smearedCentralValues.push_back(weightFromU + centralValue[col]);
      }
  return smearedCentralValues;
  } // WeightCalc::MultiGaussianSmearing() - Gives a single set of smeared parameters


  /////////
  //
  /*
    static std::vector<double>              MultiGaussianSmearing(
    std::vector<double> const& centralValue,
    TMatrixD* const& LowerTriangleCovarianceMatrix,
    bool isDecomposed,
    std::vector<double> rand);
  */
  ////////
  std::vector<double> WeightCalc::MultiGaussianSmearing(std::vector<double> const& centralValue, TMatrixD* const& LowerTriangleCovarianceMatrix, bool isDecomposed, std::vector< double > rand)
  {
    
    //compute the smeared central values
    std::vector<double> smearedCentralValues;
    
    // This guards against accidentally 
    if(!isDecomposed){
      throw art::Exception(art::errors::StdException)
	<< "Must supply the decomposed lower triangular covariance matrix.";
      return smearedCentralValues;      
    } 

     
    for(unsigned int col = 0; col < centralValue.size(); ++col)
      {
	//find the weight of each col of the upper triangular cov. matrix
	double weightFromU = 0.;
	
	for(unsigned int row = 0; row < col+1; ++row)
	  {
	    weightFromU += LowerTriangleCovarianceMatrix[0][row][col]*rand[row];
 	  }

	//multiply this weight by each col of the central values to obtain
	//the gaussian smeared and constrained central values
	// smearedCentralValues.push_back(weightFromU * centralValue[col]);
	smearedCentralValues.push_back(weightFromU + centralValue[col]);
      }
  return smearedCentralValues;
  } // WeightCalc::MultiGaussianSmearing() - Gives a single set of smeared parameters



} // namespace evwgh
 
}
