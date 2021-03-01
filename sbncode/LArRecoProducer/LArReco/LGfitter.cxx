#include "LGfitter.h"

namespace LGfitter {
// Numeric constants
constexpr Double_t mInvSq2Pi = 0.3989422804014; // (2 pi)^(-1/2)
constexpr Double_t mMPShift = -0.22278298;      // Landau maximum location

TF1* LGfitter::Fit(TH1* h, Double_t lbound, Double_t rbound, Double_t* fitparams, Double_t* fiterrors, Double_t* covmat) const
{

  // Fit the histogram h to a Landau function and pass back the parameters
  //    and their errors.
  //  Use TMath::Landau, but correct for MPV offset
  //  Initializes the fit parameters to reasonable values
  //  Pass back the fit parameters and their errors
  //  Return the fit function

  //  Note: to get this to work correctly, you need to tune the
  //    two constants "control parameters" in langaufun
  //    to match your application (mSC to match the gaus width and
  //    mNP to accomodate the histogram binning)!!

  //gStyle->SetOptFit(12);

  //  Fit histogram to Landau/Gaussian conv
  Char_t FunName[100];
  sprintf(FunName, "Fitfcn_%s", h->GetName());

  // We cannot pass langaufun directly as it is a non static member function
  // so we have to wrap in in this lambda function which both keeps ROOT happy and
  // allows us to pass the configurable parameters of NP and SC that is unique to each
  // instance of LGfitter, in turn allowing for this to be more easily controlled than
  // if it were a constexpr
  auto langaufunfun = [this](Double_t* x, Double_t* par) { return this->langaufun(x, par, this->mNP, this->mSC); };

  TF1* ffit = new TF1(FunName, langaufunfun, lbound, rbound, 4);
  ffit->SetParameters(fitparams[0], fitparams[1], fitparams[2], fitparams[3]);
  if (fiterrors) {
    std::cout << "LW error: " << fiterrors[0] << std::endl;
    ffit->SetParError(0, fiterrors[0]);
    if (fiterrors[0] == 0)
      ffit->FixParameter(0, fitparams[0]);
    ffit->SetParError(1, fiterrors[1]);
    ffit->SetParError(2, fiterrors[2]);
    ffit->SetParError(3, fiterrors[3]);
  }
  ffit->SetParLimits(0, 0.001, h->GetRMS());
  ffit->SetParLimits(1, lbound, rbound);
  ffit->SetParLimits(2, h->GetMaximum() / 2.f, h->GetMaximum() * 2.f);
  ffit->SetParLimits(3, 0.001, h->GetRMS());

  ffit->SetParNames("Width", "MP", "Amp", "Sigma");
  // If the bins are large w.r.t. to the rising slope, you may
  //     need to use the I option when fitting.  i.e. "IMLEVR"
  TFitResultPtr r = h->Fit(FunName, mFitOpt.c_str());

  // Check fit status
  fitparams[0] = -1000.0;
  fitparams[1] = -1000.0;
  fitparams[2] = -1000.0;
  fitparams[3] = -1000.0;
  if (fiterrors) {
    fiterrors[0] = -1000.0;
    fiterrors[1] = -1000.0;
    fiterrors[2] = -1000.0;
    fiterrors[3] = -1000.0;
  }
  if (covmat) {
    covmat[0] = -1000.0;
    covmat[1] = -1000.0;
    covmat[2] = -1000.0;
    covmat[3] = -1000.0;
  }
  //  if (ii<20) return(ffit);
  if ((Int_t)r == 0) { //successful fit
    // Get Fit Parameters, their errors and cov matrix
    fitparams[0] = h->GetFunction(FunName)->GetParameter(0);
    fitparams[1] = h->GetFunction(FunName)->GetParameter(1);
    fitparams[2] = h->GetFunction(FunName)->GetParameter(2);
    fitparams[3] = h->GetFunction(FunName)->GetParameter(3);
    if (fiterrors) {
      fiterrors[0] = h->GetFunction(FunName)->GetParError(0);
      fiterrors[1] = h->GetFunction(FunName)->GetParError(1);
      fiterrors[2] = h->GetFunction(FunName)->GetParError(2);
      fiterrors[3] = h->GetFunction(FunName)->GetParError(3);
    }
    TVirtualFitter* fitter = TVirtualFitter::GetFitter();
    if (fitter && covmat) {
      TMatrixD matrix(4, 4, fitter->GetCovarianceMatrix());
      covmat[0] = fitter->GetCovarianceMatrixElement(0, 1);
      covmat[1] = fitter->GetCovarianceMatrixElement(0, 2);
      covmat[2] = fitter->GetCovarianceMatrixElement(1, 2);
      covmat[3] = fitter->GetCovarianceMatrixElement(0, 3);
      // std::cout << "cov int " << covmat[3] << std::endl;
      //missing covariance terms here !
    } else {
      std::cout << "Fitter not found" << std::endl;
    }
  } else {
    std::cout << "Fit status: " << (Int_t)r << std::endl;
  }

  return (ffit);
}

//Double_t LandFun(Double_t* x, Double_t* par)
//{

//  //Fit parameters:
//  //par[0]=Width (scale) parameter of Landau density (sigma)
//  //par[1]=Most Probable (MP, location) parameter of Landau density
//  //par[2]=Landau Amplitude
//  //
//  //In the Landau distribution (represented by the CERNLIB approximation),
//  //the maximum is located at x=-0.22278298 with the location parameter=0.
//  //This shift is corrected within this function, so that the actual
//  //maximum is identical to the MP parameter

//  const Double_t mpc(par[1] - mMPShift * par[0]);

//  const Double_t temp(par[2] * (TMath::Landau(x[0], mpc, par[0])));

//  return (temp);
//}

Double_t LGfitter::langaufun(Double_t* x, Double_t* par, const Double_t NP, const Double_t SC) const
{

  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation),
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.

  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow, xupp;
  Double_t step;
  Double_t i;

  // MP shift correction
  mpc = par[1] - mMPShift * par[0];

  // Range of convolution integral
  xlow = x[0] - SC * par[3];
  xupp = x[0] + SC * par[3];

  step = (xupp - xlow) / NP;

  if (par[0] == 0)
    sum = 0;
  else {
    // Convolution integral of Landau and Gaussian by sum
    for (i = 1.0; i <= NP / 2; i++) {
      xx = xlow + (i - .5) * step;
      fland = TMath::Landau(xx, mpc, par[0]) / par[0];
      sum += fland * TMath::Gaus(x[0], xx, par[3]);

      xx = xupp - (i - .5) * step;
      fland = TMath::Landau(xx, mpc, par[0]) / par[0];
      sum += fland * TMath::Gaus(x[0], xx, par[3]);
    }
  }

  return (par[2] * step * sum * mInvSq2Pi / par[3]);
}
}
