#pragma once

#include "CAFAna/Experiment/IExperiment.h"

#include "TMath.h"

#include "TRandom3.h"

namespace ana{class IOscCalculatorAdjustable;}

namespace ana
{
  // http://www.nu-fit.org/?q=node/139
  const double kNuFitDmsq21CV = 7.50e-5;
  const double kNuFitTh12CV = 33.56 * TMath::Pi()/180;

  // Have to adjust for nu-fit's weird convention in NH
  const double kNuFitDmsq32CVNH = +2.524e-3 - kNuFitDmsq21CV;
  const double kNuFitTh23CVNH = 41.6 * TMath::Pi()/180;
  const double kNuFitTh13CVNH = 8.46 * TMath::Pi()/180;
  const double kNuFitdCPCVNH = 261 * TMath::Pi()/180;

  const double kNuFitDmsq32CVIH = -2.514e-3;
  const double kNuFitTh23CVIH = 50.0 * TMath::Pi()/180;
  const double kNuFitTh13CVIH = 8.49 * TMath::Pi()/180;
  const double kNuFitdCPCVIH = 277 * TMath::Pi()/180;

  // Based on 1/3 of the 3sigma error
  const double kNuFitDmsq21Err = ((8.90-7.50)/3)*1e-5;
  const double kNuFitTh12Err = ((35.99-33.56)/3) * TMath::Pi()/180;

  const double kNuFitDmsq32ErrNH = ((2.643-2.524)/3)*1e-3;
  const double kNuFitTh23ErrNH = ((52.8-41.6)/3) * TMath::Pi()/180;
  const double kNuFitTh13ErrNH = ((8.90-8.46)/3) * TMath::Pi()/180;

  const double kNuFitDmsq32ErrIH = ((2.635-2.514)/3)*1e-3;
  const double kNuFitTh23ErrIH = ((53.1-50.0)/3) * TMath::Pi()/180;
  const double kNuFitTh13ErrIH = ((8.93-8.49)/3) * TMath::Pi()/180;


  // hie = +/-1
  osc::IOscCalculatorAdjustable* NuFitOscCalcCDR(int hie);

  osc::IOscCalculatorAdjustable* NuFitOscCalcCDRPlusOneSigma(int hie);

  // Add in a throw for toys
  osc::IOscCalculatorAdjustable* ThrownNuFitOscCalcCDR(int hie);

  class NuFitPenalizerCDR: public IExperiment
  {
  public:
    double ChiSq(osc::IOscCalculatorAdjustable* calc,
                 const SystShifts& syst = SystShifts::Nominal()) const override;
  };

  class Penalizer_GlbLikeCDR: public IExperiment
  {
  public:
    Penalizer_GlbLikeCDR(osc::IOscCalculatorAdjustable* cvcalc, int hietrue, bool weakOnly=false);

    double Dmsq21CV() const {return fDmsq21;}
    double Th12CV() const {return fTh12;}
    double Dmsq32CV() const {return fDmsq32;}
    double Th23CV() const {return fTh23;}
    double Th13CV() const {return fTh13;}

    double ChiSq(osc::IOscCalculatorAdjustable* calc,
                 const SystShifts& syst = SystShifts::Nominal()) const override;

  protected:
    double fDmsq21;
    double fTh12;
    double fDmsq32;
    double fTh23;
    double fTh13;
    double fRho;

    double fDmsq21Err;
    double fTh12Err;
    double fDmsq32Err;
    double fTh23Err;
    double fTh13Err;
    double fRhoErr;

  private:
    // Okay, I'm bad at naming things. This is a flag to apply a penalty to all parameters
    // or only those that are weakly constained in DUNE (12 sector, rho)
    bool fWeakOnly;
  };

}
