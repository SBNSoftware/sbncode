#ifndef _MultiVariateRNG_h_
#define _MultiVariateRNG_h_

#include "TMatrixDSym.h"
#include "TRandom2.h"

#include <vector>
#include <iostream>

#include "cetlib_except/exception.h"

const int MAX_CALLS = 200000;

class MultiVariateRNG {

 public:

  MultiVariateRNG(unsigned int seed,TMatrixDSym cov,std::vector<double> cv = std::vector<double>());
  ~MultiVariateRNG(); 

  std::vector<double> GetParameterSet() const;

 private:

  unsigned int fSeed;
  TMatrixDSym fCov;
  TMatrixDSym fCovInv;
  double fDeterminant;
  std::vector<double> fCV;

  TRandom2 *R;

};

#endif
