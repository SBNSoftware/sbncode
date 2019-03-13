#include "CAFAna/Analysis/CalcsNuFit.h"

#include "Utilities/func/MathUtil.h"

#include "OscLib/func/OscCalculatorPMNSOpt.h"

#include "CAFAna/Vars/FitVars.h"

namespace ana
{
  //----------------------------------------------------------------------
  osc::IOscCalculatorAdjustable* NuFitOscCalc(int hie)
  {
    assert(hie == +1 || hie == -1);

    osc::IOscCalculatorAdjustable* ret = new osc::OscCalculatorPMNSOpt;
    ret->SetL(kBaseline);
    ret->SetRho(kEarthDensity);

    ret->SetDmsq21(kNuFitDmsq21CV);
    ret->SetTh12(kNuFitTh12CV);

    if(hie > 0){
      ret->SetDmsq32(kNuFitDmsq32CVNH);
      ret->SetTh23(kNuFitTh23CVNH);
      ret->SetTh13(kNuFitTh13CVNH);
      ret->SetdCP(kNuFitdCPCVNH);
    }
    else{
      ret->SetDmsq32(kNuFitDmsq32CVIH);
      ret->SetTh23(kNuFitTh23CVIH);
      ret->SetTh13(kNuFitTh13CVIH);
      ret->SetdCP(kNuFitdCPCVIH);
    }

    return ret;
  }

  bool HasVar(std::vector<const IFitVar*> oscVars, std::string name){
    for(auto *s :oscVars ) if(s->ShortName() == name) return true;
    return false;
  }

  //----------------------------------------------------------------------
  osc::IOscCalculatorAdjustable* ThrownNuFitOscCalc(int hie, std::vector<const IFitVar*> oscVars)
  {
    assert(hie == +1 || hie == -1);

    osc::IOscCalculatorAdjustable* ret = NuFitOscCalc(hie);//new osc::OscCalculatorPMNSOpt;

    // Throw 12 and rho within errors
    if (HasVar(oscVars, kFitRho.ShortName()))
      ret->SetRho(kEarthDensity*(1+0.02*gRandom->Gaus()));

    if (HasVar(oscVars, kFitDmSq21.ShortName()))
      ret->SetDmsq21(kNuFitDmsq21CV*(1+kNuFitDmsq21Err*gRandom->Gaus()));

    if (HasVar(oscVars, kFitSinSq2Theta12.ShortName()))
      ret->SetTh12(kNuFitTh12CV*(1+kNuFitTh12Err*gRandom->Gaus()));

    // Uniform throws within +/-3 sigma
    if(hie > 0){
      if (HasVar(oscVars, kFitDmSq32Scaled.ShortName()))
	ret->SetDmsq32(gRandom->Uniform(kNuFitDmsq32CVNH-3*kNuFitDmsq32ErrNH,
					kNuFitDmsq32CVNH+3*kNuFitDmsq32ErrNH));

      if (HasVar(oscVars, kFitSinSqTheta23.ShortName()))
        ret->SetTh23(gRandom->Uniform(kNuFitTh23CVNH-3*kNuFitTh23ErrNH,
                                      kNuFitTh23CVNH+3*kNuFitTh23ErrNH));

      if (HasVar(oscVars, kFitTheta13.ShortName()))
        ret->SetTh13(gRandom->Uniform(kNuFitTh13CVNH-3*kNuFitTh13ErrNH,
                                      kNuFitTh13CVNH+3*kNuFitTh13ErrNH));

      if (HasVar(oscVars, kFitDeltaInPiUnits.ShortName()))
        ret->SetdCP(gRandom->Uniform(-1*TMath::Pi(), TMath::Pi()));

    } else {
      if (HasVar(oscVars, kFitDmSq32Scaled.ShortName()))
	ret->SetDmsq32(gRandom->Uniform(kNuFitDmsq32CVIH-3*kNuFitDmsq32ErrIH,
					kNuFitDmsq32CVIH+3*kNuFitDmsq32ErrIH));

      if (HasVar(oscVars, kFitSinSqTheta23.ShortName()))
        ret->SetTh23(gRandom->Uniform(kNuFitTh23CVIH-3*kNuFitTh23ErrIH,
                                      kNuFitTh23CVIH+3*kNuFitTh23ErrIH));

      if (HasVar(oscVars, kFitTheta13.ShortName()))
        ret->SetTh13(gRandom->Uniform(kNuFitTh13CVIH-3*kNuFitTh13ErrIH,
                                      kNuFitTh13CVIH+3*kNuFitTh13ErrIH));

      if (HasVar(oscVars, kFitDeltaInPiUnits.ShortName()))
	ret->SetdCP(gRandom->Uniform(-1*TMath::Pi(), TMath::Pi()));
    }
    return ret;
  }

  //----------------------------------------------------------------------
  osc::IOscCalculatorAdjustable* NuFitOscCalcPlusOneSigma(int hie)
  {
    assert(hie == +1 || hie == -1);

    osc::IOscCalculatorAdjustable* ret = new osc::OscCalculatorPMNSOpt;
    ret->SetL(kBaseline);
    ret->SetRho(kEarthDensity);

    ret->SetDmsq21(kNuFitDmsq21CV + kNuFitDmsq21Err);
    ret->SetTh12(kNuFitTh12CV + kNuFitTh12Err);

    if(hie > 0){
      ret->SetDmsq32(kNuFitDmsq32CVNH + kNuFitDmsq32ErrNH);
      ret->SetTh23(kNuFitTh23CVNH + kNuFitTh23ErrNH);
      ret->SetTh13(kNuFitTh13CVNH + kNuFitTh13ErrNH);
    }
    else{
      ret->SetDmsq32(kNuFitDmsq32CVIH + kNuFitDmsq32ErrIH);
      ret->SetTh23(kNuFitTh23CVIH + kNuFitTh23ErrIH);
      ret->SetTh13(kNuFitTh13CVIH + kNuFitTh13ErrIH);
    }

    ret->SetdCP(0); // a little ambiguous in the instructions

    return ret;
  }

  //----------------------------------------------------------------------
  double NuFitPenalizer::ChiSq(osc::IOscCalculatorAdjustable* calc,
                               const SystShifts& /*syst*/) const
  {
    double ret =
      util::sqr((calc->GetDmsq21() - kNuFitDmsq21CV)/kNuFitDmsq21Err) +
      util::sqr((calc->GetTh12() - kNuFitTh12CV)/kNuFitTh12Err);

    if(calc->GetDmsq32() > 0){
      ret +=
        util::sqr((calc->GetDmsq32() - kNuFitDmsq32CVNH)/kNuFitDmsq32ErrNH) +
        util::sqr((calc->GetTh23() - kNuFitTh23CVNH)/kNuFitTh23ErrNH) +
        util::sqr((calc->GetTh13() - kNuFitTh13CVNH)/kNuFitTh13ErrNH);
    }
    else{
      ret +=
        util::sqr((calc->GetDmsq32() - kNuFitDmsq32CVIH)/kNuFitDmsq32ErrIH) +
        util::sqr((calc->GetTh23() - kNuFitTh23CVIH)/kNuFitTh23ErrIH) +
        util::sqr((calc->GetTh13() - kNuFitTh13CVIH)/kNuFitTh13ErrIH);
    }

    // No term in delta

    return ret;
  }

  //----------------------------------------------------------------------
  Penalizer_GlbLike::Penalizer_GlbLike(osc::IOscCalculatorAdjustable* cvcalc, int hietrue, bool weakOnly) : fWeakOnly(weakOnly) {

    fDmsq21 = cvcalc->GetDmsq21();
    fTh12 = cvcalc->GetTh12();
    fDmsq32 = cvcalc->GetDmsq32();
    fTh23 = cvcalc->GetTh23();
    fTh13 = cvcalc->GetTh13();
    fRho = cvcalc->GetRho();

    //Set the errors by hand for now.
    //Fractional errors in GLoBES convention
    //NH: 0.023, 0.018, 0.058, 0.0, 0.024, 0.016
    //IH: 0.023, 0.018, 0.048, 0.0, 0.024, 0.016

    fDmsq21Err = kNuFitDmsq21Err;
    fTh12Err = kNuFitTh12Err;
    fDmsq32Err = (hietrue > 0) ? kNuFitDmsq32ErrNH : kNuFitDmsq32ErrIH;
    fTh13Err = (hietrue > 0) ? kNuFitTh13ErrNH : kNuFitTh13ErrIH;
    fTh23Err = (hietrue > 0) ? kNuFitTh23ErrNH : kNuFitTh23ErrIH;

    fRhoErr = 0.02*fRho;

  }

  double Penalizer_GlbLike::ChiSq(osc::IOscCalculatorAdjustable* calc,
				  const SystShifts& /*syst*/) const {

  //Usage: calc is the fit parameters as above
  //Starting fit parameters and errors are set in constructor - this is equivalent to SetCentralValues in globes

    double ret =
      util::sqr((calc->GetDmsq21() - fDmsq21)/fDmsq21Err) +
      util::sqr((calc->GetTh12() - fTh12)/fTh12Err) +
      util::sqr((calc->GetRho() - fRho)/fRhoErr);

    // if fWeakOnly is set, only apply a constraint to the parameter we can only weakly constrain in DUNE
    if (!fWeakOnly)
      ret +=
	util::sqr((calc->GetDmsq32() - fDmsq32)/fDmsq32Err) +
	util::sqr((calc->GetTh23() - fTh23)/fTh23Err) +
	util::sqr((calc->GetTh13() - fTh13)/fTh13Err);

    // No term in delta

    return ret;
  }


}
