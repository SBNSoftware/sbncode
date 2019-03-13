#include "CAFAna/Analysis/CalcsNuFit_cdr.h"

#include "Utilities/func/MathUtil.h"

#include "OscLib/func/OscCalculatorPMNSOpt.h"

namespace ana
{
  //----------------------------------------------------------------------
  osc::IOscCalculatorAdjustable* NuFitOscCalcCDR(int hie)
  {
    assert(hie == +1 || hie == -1);

    osc::IOscCalculatorAdjustable* ret = new osc::OscCalculatorPMNSOpt;
    ret->SetL(1300);
    ret->SetRho(2.95674); // g/cm^3. Dan Cherdack's doc "used in GLOBES"

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

  //----------------------------------------------------------------------
  osc::IOscCalculatorAdjustable* ThrownNuFitOscCalcCDR(int hie)
  {
    assert(hie == +1 || hie == -1);

    osc::IOscCalculatorAdjustable* ret = new osc::OscCalculatorPMNSOpt;
    ret->SetL(1300);

    // Throw 12 and rho within errors
    ret->SetRho(2.95674*(1+0.02*gRandom->Gaus()));
    ret->SetDmsq21(kNuFitDmsq21CV*(1+kNuFitDmsq21Err*gRandom->Gaus()));
    ret->SetTh12(kNuFitTh12CV*(1+kNuFitTh12Err*gRandom->Gaus()));

    // Uniform throws within +/-3 sigma
    if(hie > 0){
      ret->SetDmsq32(gRandom->Uniform(kNuFitDmsq32CVNH-3*kNuFitDmsq32ErrNH, 
				      kNuFitDmsq32CVNH+3*kNuFitDmsq32ErrNH));
      ret->SetTh23(gRandom->Uniform(kNuFitTh23CVNH-3*kNuFitTh23ErrNH,
				    kNuFitTh23CVNH+3*kNuFitTh23ErrNH));
      ret->SetTh13(gRandom->Uniform(kNuFitTh13CVNH-3*kNuFitTh13ErrNH,
				    kNuFitTh13CVNH+3*kNuFitTh13ErrNH));
    } else {
      ret->SetDmsq32(gRandom->Uniform(kNuFitDmsq32CVIH-3*kNuFitDmsq32ErrIH,
                                      kNuFitDmsq32CVIH+3*kNuFitDmsq32ErrIH));
      ret->SetTh23(gRandom->Uniform(kNuFitTh23CVIH-3*kNuFitTh23ErrIH,
                                    kNuFitTh23CVIH+3*kNuFitTh23ErrIH));
      ret->SetTh13(gRandom->Uniform(kNuFitTh13CVIH-3*kNuFitTh13ErrIH,
                                    kNuFitTh13CVIH+3*kNuFitTh13ErrIH));
    }
    ret->SetdCP(gRandom->Uniform(-1*TMath::Pi(), TMath::Pi()));
		  
    return ret;
  }


  //----------------------------------------------------------------------
  osc::IOscCalculatorAdjustable* NuFitOscCalcCDRPlusOneSigma(int hie)
  {
    assert(hie == +1 || hie == -1);

    osc::IOscCalculatorAdjustable* ret = new osc::OscCalculatorPMNSOpt;
    ret->SetL(1300);
    ret->SetRho(2.95674); // g/cm^3. Dan Cherdack's doc "used in GLOBES"

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
  double NuFitPenalizerCDR::ChiSq(osc::IOscCalculatorAdjustable* calc,
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
  Penalizer_GlbLikeCDR::Penalizer_GlbLikeCDR(osc::IOscCalculatorAdjustable* cvcalc, int hietrue, bool weakOnly) : fWeakOnly(weakOnly) {

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

    fDmsq21Err = 0.024*fDmsq21;
    fTh12Err = 0.023*fTh12;
    fDmsq32Err = 0.016*fDmsq32;
    fTh13Err = 0.018*fTh13;
    fTh23Err = (hietrue > 0) ? 0.058*fTh23 : 0.048*fTh23;

    fRhoErr = 0.02*fRho;
    
  }
    
  double Penalizer_GlbLikeCDR::ChiSq(osc::IOscCalculatorAdjustable* calc,
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
