#pragma once

#include "CAFAna/Experiment/IExperiment.h"
#include "CAFAna/Core/IFitVar.h"

#include "TMath.h"

#include "TRandom3.h"

namespace ana{class IOscCalculatorAdjustable;}

namespace ana
{
  // http://www.nu-fit.org/?q=node/177
  // NuFit November 2018
  const double kNuFitDmsq21CV = 7.39e-5;
  const double kNuFitTh12CV = 33.82 * TMath::Pi()/180;

  // Have to adjust for nu-fit's weird convention in NH
  const double kNuFitDmsq32CVNH = +2.525e-3 - kNuFitDmsq21CV;
  const double kNuFitTh23CVNH = 49.6 * TMath::Pi()/180;
  const double kNuFitTh13CVNH = 8.61 * TMath::Pi()/180;
  const double kNuFitdCPCVNH = 215 * TMath::Pi()/180;

  const double kNuFitDmsq32CVIH = -2.512e-3;
  const double kNuFitTh23CVIH = 49.8 * TMath::Pi()/180;
  const double kNuFitTh13CVIH = 8.65 * TMath::Pi()/180;
  const double kNuFitdCPCVIH = 284 * TMath::Pi()/180;

  // Based on 1/6 of the +/- 3sigma error
  const double kNuFitDmsq21Err = ((8.01-6.79)/6)*1e-5;
  const double kNuFitTh12Err = ((36.27-31.61)/6) * TMath::Pi()/180;

  const double kNuFitDmsq32ErrNH = ((2.625-2.427)/6)*1e-3;
  const double kNuFitTh23ErrNH = ((52.4-40.3)/6) * TMath::Pi()/180;
  const double kNuFitTh13ErrNH = ((8.99-8.22)/6) * TMath::Pi()/180;

  const double kNuFitDmsq32ErrIH = ((2.611-2.412)/6)*1e-3;
  const double kNuFitTh23ErrIH = ((52.5-40.6)/6) * TMath::Pi()/180;
  const double kNuFitTh13ErrIH = ((9.03-8.27)/6) * TMath::Pi()/180;

  //https://arxiv.org/pdf/1707.02322.pdf
  const double kBaseline = 1284.9;     // km
  const double kEarthDensity = 2.848;  // g/cm^3

  // hie = +/-1
  osc::IOscCalculatorAdjustable* NuFitOscCalc(int hie);

  osc::IOscCalculatorAdjustable* NuFitOscCalcPlusOneSigma(int hie);

  // Add in a throw for toys
  osc::IOscCalculatorAdjustable* ThrownNuFitOscCalc(int hie, std::vector<const IFitVar*> oscVars);

  bool HasVar(std::vector<const IFitVar*> oscVars, std::string name);


  class NuFitPenalizer: public IExperiment
  {
  public:
    double ChiSq(osc::IOscCalculatorAdjustable* calc,
                 const SystShifts& syst = SystShifts::Nominal()) const override;

    void Derivative(osc::IOscCalculator*,
                    const SystShifts&,
                    std::unordered_map<const ISyst*, double>&) const override
    {
      // Nothing to do, since we only depend on osc params and the derivatives
      // with the systs are all zero. But need to implement, because the
      // default implementation indicates that we are unable to calculate the
      // necessary derivatives.
    }
  };

  class Penalizer_GlbLike: public IExperiment
  {
  public:
    Penalizer_GlbLike(osc::IOscCalculatorAdjustable* cvcalc, int hietrue, bool weakOnly=false);

    double Dmsq21CV() const {return fDmsq21;}
    double Th12CV() const {return fTh12;}
    double Dmsq32CV() const {return fDmsq32;}
    double Th23CV() const {return fTh23;}
    double Th13CV() const {return fTh13;}
    double RhoCV() const {return fRho;}

    double Dmsq21Err() const {return fDmsq21Err;}
    double Th12Err() const {return fTh12Err;}
    double Dmsq32Err() const {return fDmsq32Err;}
    double Th23Err() const {return fTh23Err;}
    double Th13Err() const {return fTh13Err;}
    double RhoErr() const {return fRhoErr;}

    double ChiSq(osc::IOscCalculatorAdjustable* calc,
                 const SystShifts& syst = SystShifts::Nominal()) const override;

    void Derivative(osc::IOscCalculator*,
                    const SystShifts&,
                    std::unordered_map<const ISyst*, double>&) const override
    {
      // See justification in NuFitPenalizer::Derivative()
    }

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
