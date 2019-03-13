#include "CAFAna/Vars/FitVarsSterile.h"

#include "OscLib/func/IOscCalculator.h"
#include "OscLib/func/OscCalculatorSterile.h"
#include "Utilities/func/MathUtil.h"

#include <cassert>
#include <cmath>

namespace ana
{
  
  //--------------------------------------------------------------------------- 
  double FitDmSq32Sterile::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    const osc::OscCalculatorSterile* sterile = osc::DowncastToSterile(osc);
    double dm221 = sterile->GetDm(2);
    double dm231 = sterile->GetDm(3);
    return dm231 - dm221;
  }
  
  //---------------------------------------------------------------------------                   
  void FitDmSq32Sterile::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc::OscCalculatorSterile* sterile = osc::DowncastToSterile(osc);
    double dm221 = sterile->GetDm(2);
    double dm231 = val + dm221;
    sterile->SetDm(3, dm231);
  }

  //---------------------------------------------------------------------------                   
  double FitDmSq41Sterile::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return osc::DowncastToSterile(osc)->GetDm(4);
  }
  
  //---------------------------------------------------------------------------                   
  void FitDmSq41Sterile::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc::DowncastToSterile(osc)->SetDm(4, val);
  }

  //--------------------------------------------------------------------------- 
  double FitDmSq43Sterile::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    const osc::OscCalculatorSterile* sterile = osc::DowncastToSterile(osc);
    double dm241 = sterile->GetDm(4);
    double dm231 = sterile->GetDm(3);
    return dm241 - dm231;
  }
  
  //---------------------------------------------------------------------------                   
  void FitDmSq43Sterile::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc::OscCalculatorSterile* sterile = osc::DowncastToSterile(osc);
    double dm231 = sterile->GetDm(3);
    double dm241 = val + dm231;
    sterile->SetDm(4, dm241);
  }
  
  //----------------------------------------------------------------------
  double FitDelta13InPiUnitsSterile::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    const osc::OscCalculatorSterile* sterile = osc::DowncastToSterile(osc);
    double ret = sterile->GetDelta(1, 3)/M_PI;
    while(ret < 0) ret += 2;
    while(ret > 2) ret -= 2;
    return ret;
  }

  //----------------------------------------------------------------------
  void FitDelta13InPiUnitsSterile::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc::DowncastToSterile(osc)->SetDelta(1, 3, M_PI*val);
  }

  //----------------------------------------------------------------------
  double FitDelta14InPiUnitsSterile::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    const osc::OscCalculatorSterile* sterile = osc::DowncastToSterile(osc);
    double ret = sterile->GetDelta(1, 4)/M_PI;
    while(ret < 0) ret += 2;
    while(ret > 2) ret -= 2;
    return ret;
  }

  //----------------------------------------------------------------------
  void FitDelta14InPiUnitsSterile::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc::DowncastToSterile(osc)->SetDelta(1, 4, M_PI*val);
  }

  //----------------------------------------------------------------------
  double FitDelta24InPiUnitsSterile::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    const osc::OscCalculatorSterile* sterile = osc::DowncastToSterile(osc);
    double ret = sterile->GetDelta(2, 4)/M_PI;
    while(ret < 0) ret += 2;
    while(ret > 2) ret -= 2;
    return ret;
  }

  //----------------------------------------------------------------------
  void FitDelta24InPiUnitsSterile::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc::DowncastToSterile(osc)->SetDelta(2, 4, M_PI*val);
  }

  //---------------------------------------------------------------------------                   
  double FitTheta13Sterile::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return osc::DowncastToSterile(osc)->GetAngle(1,3);
  }

  //----------------------------------------------------------------------
  void FitTheta13Sterile::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc::DowncastToSterile(osc)->SetAngle(1, 3, Clamp(val));
  }

  //----------------------------------------------------------------------
  double FitSinSqTheta13Sterile::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return util::sqr(sin(osc::DowncastToSterile(osc)->GetAngle(1,3)));
  }

  //----------------------------------------------------------------------
  void FitSinSqTheta13Sterile::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc::DowncastToSterile(osc)->SetAngle(1, 3, asin(sqrt(Clamp(val))));
  }

  //---------------------------------------------------------------------------                   
  double FitTheta23Sterile::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return osc::DowncastToSterile(osc)->GetAngle(2,3);
  }

  //----------------------------------------------------------------------
  void FitTheta23Sterile::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc::DowncastToSterile(osc)->SetAngle(2, 3, Clamp(val));
  }

  //----------------------------------------------------------------------
  double FitSinSqTheta23Sterile::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return util::sqr(sin(osc::DowncastToSterile(osc)->GetAngle(2,3)));
  }

  //----------------------------------------------------------------------
  void FitSinSqTheta23Sterile::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc::DowncastToSterile(osc)->SetAngle(2, 3, asin(sqrt(Clamp(val))));
  }

  //---------------------------------------------------------------------------                   
  double FitTheta14Sterile::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return osc::DowncastToSterile(osc)->GetAngle(1,4);
  }

  //----------------------------------------------------------------------
  void FitTheta14Sterile::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc::DowncastToSterile(osc)->SetAngle(1, 4, Clamp(val));
  }

  //----------------------------------------------------------------------
  double FitSinSqTheta14Sterile::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return util::sqr(sin(osc::DowncastToSterile(osc)->GetAngle(1,4)));
  }

  //----------------------------------------------------------------------
  void FitSinSqTheta14Sterile::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc::DowncastToSterile(osc)->SetAngle(1, 4, asin(sqrt(Clamp(val))));
  }

  //----------------------------------------------------------------------
  double FitSinSq2Theta14Sterile::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return util::sqr(sin(2*osc::DowncastToSterile(osc)->GetAngle(1,4)));
  }

  //----------------------------------------------------------------------
  void FitSinSq2Theta14Sterile::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc::DowncastToSterile(osc)->SetAngle(1,4, asin(sqrt(Clamp(val)))/2);
  }

  //---------------------------------------------------------------------------                   
  double FitTheta24Sterile::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return osc::DowncastToSterile(osc)->GetAngle(2,4);
  }

  //----------------------------------------------------------------------
  void FitTheta24Sterile::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc::DowncastToSterile(osc)->SetAngle(2, 4, Clamp(val));
  }

  //----------------------------------------------------------------------
  double FitSinSqTheta24Sterile::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return util::sqr(sin(osc::DowncastToSterile(osc)->GetAngle(2,4)));
  }

  //----------------------------------------------------------------------
  void FitSinSqTheta24Sterile::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc::DowncastToSterile(osc)->SetAngle(2, 4, asin(sqrt(Clamp(val))));
  }

  //----------------------------------------------------------------------
  double FitSinSq2Theta24Sterile::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return util::sqr(sin(2*osc::DowncastToSterile(osc)->GetAngle(2,4)));
  }

  //----------------------------------------------------------------------
  void FitSinSq2Theta24Sterile::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc::DowncastToSterile(osc)->SetAngle(2,4, asin(sqrt(Clamp(val)))/2);
  }

  //---------------------------------------------------------------------------                   
  double FitTheta34Sterile::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return osc::DowncastToSterile(osc)->GetAngle(3,4);
  }

  //----------------------------------------------------------------------
  void FitTheta34Sterile::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc::DowncastToSterile(osc)->SetAngle(3, 4, Clamp(val));
  }

  //----------------------------------------------------------------------
  double FitSinSqTheta34Sterile::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return util::sqr(sin(osc::DowncastToSterile(osc)->GetAngle(3,4)));
  }

  //----------------------------------------------------------------------
  void FitSinSqTheta34Sterile::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc::DowncastToSterile(osc)->SetAngle(3, 4, asin(sqrt(Clamp(val))));
  }

  //----------------------------------------------------------------------
  double FitSinSq2Theta34Sterile::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return util::sqr(sin(2*osc::DowncastToSterile(osc)->GetAngle(3,4)));
  }

  //----------------------------------------------------------------------
  void FitSinSq2Theta34Sterile::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc::DowncastToSterile(osc)->SetAngle(3,4, asin(sqrt(Clamp(val)))/2);
  }

  //---------------------------------------------------------------------------                   
  double FitTheta13InDegreesSterile::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return TMath::RadToDeg()*osc::DowncastToSterile(osc)->GetAngle(1,3);
  }

  //----------------------------------------------------------------------
  void FitTheta13InDegreesSterile::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc::DowncastToSterile(osc)->SetAngle(1, 3, TMath::DegToRad()*Clamp(val));
  }

  //---------------------------------------------------------------------------                   
  double FitTheta23InDegreesSterile::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return TMath::RadToDeg()*osc::DowncastToSterile(osc)->GetAngle(2,3);
  }

  //----------------------------------------------------------------------
  void FitTheta23InDegreesSterile::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc::DowncastToSterile(osc)->SetAngle(2, 3, TMath::DegToRad()*Clamp(val));
  }

  //---------------------------------------------------------------------------                   
  double FitTheta14InDegreesSterile::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return TMath::RadToDeg()*osc::DowncastToSterile(osc)->GetAngle(1,4);
  }

  //----------------------------------------------------------------------
  void FitTheta14InDegreesSterile::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc::DowncastToSterile(osc)->SetAngle(1, 4, TMath::DegToRad()*Clamp(val));
  }

  //---------------------------------------------------------------------------                   
  double FitTheta24InDegreesSterile::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return TMath::RadToDeg()*osc::DowncastToSterile(osc)->GetAngle(2,4);
  }

  //----------------------------------------------------------------------
  void FitTheta24InDegreesSterile::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc::DowncastToSterile(osc)->SetAngle(2, 4, TMath::DegToRad()*Clamp(val));
  }

  //---------------------------------------------------------------------------                   
  double FitTheta34InDegreesSterile::GetValue(const osc::IOscCalculatorAdjustable* osc) const
  {
    return TMath::RadToDeg()*osc::DowncastToSterile(osc)->GetAngle(3,4);
  }

  //----------------------------------------------------------------------
  void FitTheta34InDegreesSterile::SetValue(osc::IOscCalculatorAdjustable* osc, double val) const
  {
    osc::DowncastToSterile(osc)->SetAngle(3, 4, TMath::DegToRad()*Clamp(val));
  }

} // namespace
