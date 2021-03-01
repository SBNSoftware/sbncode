#ifndef LGFITTER_H_SEEN
#define LGFITTER_H_SEEN

//C++ Includes
#include <iostream>
#include <string>

//Root Includes
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TH1.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMinuit.h"
#include "TROOT.h"
#include "TVirtualFitter.h"

// Code for fitting a Landau-Guassian convolution to a histgram.
// Based upon the langaus.C ROOT tutorial
// Modified by M. Stancari and E. Tyley

namespace LGfitter {
class LGfitter {

  public:
  LGfitter(const Double_t np = 10000.0, const Double_t sc = 8.0, const int verbose = 0)
      : mNP(np)
      , mSC(sc)
      , mVerbose(verbose)
      , mFitOpt(mVerbose == 0 ? "LQRE" : mVerbose == 1 ? "LRE" : "LVRE")
  {
  }

  TF1* Fit(TH1* h, Double_t lbound, Double_t rbound, Double_t* fitparams, Double_t* fiterrors, Double_t* covmat) const;
  // Double_t LandFun(Double_t* x, Double_t* par) const;
  Double_t langaufun(Double_t* x, Double_t* par, const Double_t NP, const Double_t SC) const;

  private:
  // Control constants - need to be tuned for each specific application
  const Double_t mNP;        // Number of convolution steps
  const Double_t mSC;        // Convolution extends to +-sc Gaussian sigmas
  const int mVerbose;        // Level of verbosity
  const std::string mFitOpt; // Fit options passed to ROOT
};
}

#endif
