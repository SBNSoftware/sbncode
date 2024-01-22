#ifndef _MultiVariateRNG_cxx_
#define _MultiVariateRNG_cxx_

#include "MultiVariateRNG.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////////

MultiVariateRNG::MultiVariateRNG(unsigned int seed,TMatrixDSym cov,std::vector<double> cv):
  fSeed(seed),
  fCov(cov.GetNcols()),
  fCovInv(cov.GetNcols())
{   

  fCov = cov;
  fDeterminant = fCov.Determinant();

  fCovInv = fCov;
  fCovInv.Invert();

  if(!cv.size()){
    std::cout << "Central value vector not provided, assuming (1,1,1,...,1)" << std::endl;
    fCV = std::vector<double>(fCov.GetNcols(),1.0);
  }

  else if(cv.size() != (size_t)fCov.GetNcols()) 
    throw std::invalid_argument("Central value vector and covariance matrix have difference sizes!");

  else fCV = cv;

  R = new TRandom2(seed);

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

MultiVariateRNG::~MultiVariateRNG(){
  delete R;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

std::vector<double> MultiVariateRNG::GetParameterSet() const {

  int call=0;

  while(call < MAX_CALLS){

    // Generate a test vector
    std::vector<double> vals;

    for(size_t i_p=0;i_p<fCV.size();i_p++)
      vals.push_back(R->Uniform(0.0,fCV.at(i_p)+3*sqrt(fCov(i_p,i_p))));

    // Calculate the height of the multivar Gaussian 
    // at this point in the parameter space
    double exp_arg=0.0;

    for(Int_t i=0;i<fCov.GetNcols();i++)
      for(Int_t j=0;j<fCov.GetNcols();j++)
	exp_arg += (vals.at(i) - fCV.at(i))*fCovInv(i,j)*(vals.at(j) - fCV.at(j));

    double denominator = sqrt( pow(2*3.1416 , fCov.GetNcols())*fDeterminant );
    double PD = exp(-0.5*exp_arg)/denominator;
    double value = R->Uniform(0,1.0/denominator);

    if(value < PD) return vals;

    call++;
  }

  throw cet::exception("MultiVariateRNG") << "Call limit of " << MAX_CALLS << " reached" << std::endl;

  return {};
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////

#endif
