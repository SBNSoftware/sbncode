//////////////////////////////////////////////////////////////////////////
//                                                                      //
// \file OscCalculator.cxx                                              //
//                                                                      //
// Class with methods for calculating all things related to oscillation //
// probabilities.                                                       //
// <rbpatter@caltech.edu>						//
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include "OscLib/OscCalculator.h"

#include <iostream>
#include <cmath>

#include "TF1.h"
#include "TMath.h"

namespace osc {

  // --------------------------------------------
  OscCalculator::OscCalculator()
  { 

    // put some sensible defaults here...

    fRho = 2.75; // g/cm^3
    fL = 810; // km
    fDmsq21 = 7.59E-5; // eV^2
    fDmsq32 = 2.43E-3; //eV^2
    fTh12 = 0.601;
    fTh13 = 0.0;
    fTh23 = 7.85398163397448279e-01; // pi/4
    fdCP = 0;

    fUpdated = false;
  }

  // --------------------------------------------
  OscCalculator::~OscCalculator()
  {
  }

  // --------------------------------------------
  IOscCalculatorAdjustable* OscCalculator::Copy() const
  {
    return new OscCalculator(*this);
  }

  // --------------------------------------------
  double OscCalculator::P(int flavBefore, int flavAfter, double E)
  {
    bool antinu = (flavBefore<0&&flavAfter<0);
    if (antinu) {
      flavBefore *= -1;
      flavAfter  *= -1;
    }
    if      (flavBefore==12&&flavAfter==12) return P_ee(E,antinu);
    else if (flavBefore==12&&flavAfter==14) return P_em(E,antinu);
    else if (flavBefore==12&&flavAfter==16) return P_et(E,antinu);
    else if (flavBefore==14&&flavAfter==12) return P_me(E,antinu);
    else if (flavBefore==14&&flavAfter==14) return P_mm(E,antinu);
    else if (flavBefore==14&&flavAfter==16) return P_mt(E,antinu);
    else if (flavBefore==16&&flavAfter==12) return P_te(E,antinu);
    else if (flavBefore==16&&flavAfter==14) return P_tm(E,antinu);
    else if (flavBefore==16&&flavAfter==16) return P_tt(E,antinu);
    else return 0;
  }

  double OscCalculator::P_ee(double E, bool antinu) { return P_internal_ee(E,antinu,0); }
  double OscCalculator::P_em(double E, bool antinu) { return P_internal_me(E,antinu,1); }
  double OscCalculator::P_et(double E, bool antinu) { return P_internal_te(E,antinu,1); }

  double OscCalculator::P_me(double E, bool antinu) { return P_internal_me(E,antinu,0); }
  double OscCalculator::P_mm(double E, bool antinu) { return 1-P_me(E,antinu)-P_mt(E,antinu); }
  double OscCalculator::P_mt(double E, bool antinu) { return P_internal_mt(E,antinu,0); }

  double OscCalculator::P_te(double E, bool antinu) { return P_internal_te(E,antinu,0); }
  double OscCalculator::P_tm(double E, bool antinu) { return P_internal_mt(E,antinu,1); }
  double OscCalculator::P_tt(double E, bool antinu) { return 1 - P_te(E,antinu) - P_tm(E,antinu); }

  // --------------------------------------------
  TF1* OscCalculator::GetTF1(int flavBefore, int flavAfter)
  {
    TF1 *theTF1 = new TF1(Form("OscCalculatorFunction_%d_%d_%p",flavBefore,flavAfter,(void*)this),
			  this,&osc::OscCalculator::P_wrapper,0,120,2,"OscCalculator","P_wrapper");
    theTF1->SetParameters(flavBefore,flavAfter);
    theTF1->SetNpx(1000);
    return theTF1;
  }

  // --------------------------------------------
  double OscCalculator::P_wrapper(double *x, double *p)
  {
    // function for use by TF1
    int flavBefore = int(p[0]);
    int flavAfter  = int(p[1]);
    bool antinu = (flavBefore<0&&flavAfter<0);
    if (antinu) {
      flavBefore *= -1;
      flavAfter  *= -1;
    } 

    double (osc::OscCalculator::*P_xx)(double,bool);
    if      (flavBefore==12&&flavAfter==12) P_xx = &osc::OscCalculator::P_ee;
    else if (flavBefore==12&&flavAfter==14) P_xx = &osc::OscCalculator::P_em;
    else if (flavBefore==12&&flavAfter==16) P_xx = &osc::OscCalculator::P_et;
    else if (flavBefore==14&&flavAfter==12) P_xx = &osc::OscCalculator::P_me;
    else if (flavBefore==14&&flavAfter==14) P_xx = &osc::OscCalculator::P_mm;
    else if (flavBefore==14&&flavAfter==16) P_xx = &osc::OscCalculator::P_mt;
    else if (flavBefore==16&&flavAfter==12) P_xx = &osc::OscCalculator::P_te;
    else if (flavBefore==16&&flavAfter==14) P_xx = &osc::OscCalculator::P_tm;
    else if (flavBefore==16&&flavAfter==16) P_xx = &osc::OscCalculator::P_tt;
    else P_xx = &osc::OscCalculator::P_null;

    return (this->*P_xx)(x[0],antinu);
  }

  // -----------------------------------------------
  void OscCalculator::UpdateBasic()
  {
    if (fUpdated) return;

    fDmsq31 = fDmsq21 + fDmsq32;
    if (fDmsq31!=0) {
      falpha =fDmsq21 / fDmsq31;
    }
    else {
      std::cerr << "OscCalculator::UpdateBasic() -- fDmsq31 should never be zero, but it is" << std::endl;
      falpha = 0;
    }
    fsin_th12 = sin(fTh12);
    fsin_th13 = sin(fTh13);
    fsin_th23 = sin(fTh23);
    fcos_th12 = cos(fTh12);
    fcos_th13 = cos(fTh13);
    fcos_th23 = cos(fTh23);
    fsin_2th12 = sin(2*fTh12);
    fsin_2th13 = sin(2*fTh13);
    fsin_2th23 = sin(2*fTh23);
    fcos_2th12 = cos(2*fTh12);
    fcos_2th13 = cos(2*fTh13);
    fcos_2th23 = cos(2*fTh23);
    fsin_sq_th12 = fsin_th12*fsin_th12;
    fsin_sq_th13 = fsin_th13*fsin_th13;
    fsin_sq_th23 = fsin_th23*fsin_th23;
    fcos_sq_th12 = fcos_th12*fcos_th12;
    fcos_sq_th13 = fcos_th13*fcos_th13;
    fcos_sq_th23 = fcos_th23*fcos_th23;
    fsin_sq_2th12 = fsin_2th12*fsin_2th12;
    fsin_sq_2th13 = fsin_2th13*fsin_2th13;
    fsin_sq_2th23 = fsin_2th23*fsin_2th23;
    fcos_sq_2th12 = fcos_2th12*fcos_2th12;
    fcos_sq_2th13 = fcos_2th13*fcos_2th13;
    fcos_sq_2th23 = fcos_2th23*fcos_2th23;

    static const double ZperA = 0.5; // e- per nucleon
    static const double G_F = 1.16637E-23; // eV^-2
    static const double hbar_c_eV_cm = 1.97326938E-5; // eV-cm

    fV = TMath::Sqrt2()*G_F*fRho*ZperA*TMath::Na()*hbar_c_eV_cm*hbar_c_eV_cm*hbar_c_eV_cm;

    fUpdated = true;
  }

  // --------------------------------------------
  void OscCalculator::UpdateEDep(double E, bool antinu, bool fliptime)
  {
    static const double hbar_c_eV_km = 1.97326938E-10; // eV-km
    static const double eVPerGeV = 1E9;

    int s = (antinu)?-1:1;
    int t = (fliptime)?-1:1;

    fA = s*2*fV*E*eVPerGeV/fDmsq31;
    fD = fDmsq31*fL/(4*E*eVPerGeV*hbar_c_eV_km);

    fdCPproxy = s*t*fdCP;
    fsin_dCPproxy = sin(fdCPproxy);
    fcos_dCPproxy = cos(fdCPproxy);

    if (falpha!=0) {
      fC12 = TMath::Sqrt(fsin_sq_2th12+(fcos_2th12 - fA/falpha)*(fcos_2th12 - fA/falpha));
    }
    else {
      std::cerr << "OscCalculator::UpdateEDep() -- falpha should never be zero, but it is" << std::endl;
      fC12 = 1;
    }
    fC13 = TMath::Sqrt(fsin_sq_2th13+(fA-fcos_2th13)*(fA-fcos_2th13));
  }

  // --------------------------------------------
  double OscCalculator::P_internal_me(double E, bool antinu, bool fliptime)
  {
    UpdateBasic();
    UpdateEDep(E,antinu,fliptime);

    double cosC13D = cos(fC13*fD);
    double sinC13D = sin(fC13*fD);
    double sin1pAD = sin((fA+1)*fD);
    double cos1pAD = cos((fA+1)*fD);
    double sinAD = sin(fA*fD);
    double sinAm1D = sin((fA-1)*fD);
    double cosdpD = cos(fdCPproxy+fD);
    double sinApam2D = sin((fA+falpha-2)*fD);
    double cosApam2D = cos((fA+falpha-2)*fD);
    double cosaC12D = cos(falpha*fC12*fD);
    double sinaC12D = sin(falpha*fC12*fD);

    // This is coming straight from the MINOS NueAna package...

    // First we calculate the terms for the alpha expansion (good to all orders in th13)

    // Leading order term 
    double p1 = fsin_sq_th23*fsin_sq_2th13*sinC13D*sinC13D/(fC13*fC13);

    // Terms that appear at order alpha
    double p2Inner =
      fD*cosC13D*(1-fA*fcos_2th13)/fC13 - 
      fA*sinC13D*(fcos_2th13-fA)/(fC13*fC13); 

    double p2 = -2*fsin_sq_th12*fsin_sq_th23*fsin_sq_2th13*sinC13D/(fC13*fC13)*p2Inner*falpha;

    double p3Inner =
      -fsin_dCPproxy*(cosC13D - cos1pAD)*fC13 
      + fcos_dCPproxy*(fC13*sin1pAD - (1-fA*fcos_2th13)*sinC13D);

    double p3 = fsin_2th12*fsin_2th23*fsin_th13*sinC13D/(fA*fC13*fC13)*p3Inner*falpha;

    //  p1 + p2 + p3 is the complete contribution for this expansion
  
    // Now for the expansion in orders of sin(th13) (good to all order alpha) 

    double pa1 = 0.0, pa2 = 0.0;
    if (fabs(falpha)>1E-10) {
      // leading order term
      pa1 = fcos_th23*fcos_th23*fsin_sq_2th12*sinaC12D*sinaC12D/(fC12*fC12);

      // the first order in s13 term
      double t1 = (fcos_2th12 - fA/falpha)/fC12 
	- falpha*fA*fC12*fsin_sq_2th12/(2*(1-falpha)*fC12*fC12);
      double t2 = -fcos_dCPproxy*(sinApam2D-sinaC12D*t1);
      double t3 = -(cosaC12D-cosApam2D)*fsin_dCPproxy;
      double denom = (1-fA-falpha+fA*falpha*fcos_th12*fcos_th12)*fC12;
      double t4 = fsin_2th12*fsin_2th23*(1-falpha)*sinaC12D/denom;

      pa2 = t4*(t3+t2)*fsin_th13;
    }
    // pa1+pa2 is the complete contribution from this expansion

    // Now we need to add the two expansions and subtract off the terms that are
    // in both (falpha^1, s13^1)

    double t1 = sinAD*cosdpD*sinAm1D/(fA*(fA-1));
    double repeated = 2*falpha*fsin_2th12*fsin_2th23*fsin_th13*t1;

    // Calculate the total probability
    double totalP = p1+p2+p3 + (pa1+pa2) - repeated;
    return totalP;
  }

  // --------------------------------------------
  double OscCalculator::P_internal_te(double E, bool antinu, bool fliptime)
  {
    UpdateBasic();
    UpdateEDep(E,antinu,fliptime);

    double cosC13D = cos(fC13*fD);
    double sinC13D = sin(fC13*fD);
    double sin1pAD = sin((fA+1)*fD);
    double cos1pAD = cos((fA+1)*fD);
    double sinAD = sin(fA*fD);
    double sinAm1D = sin((fA-1)*fD);
    double cosdpD = cos(fdCPproxy+fD);
    double sinApam2D = sin((fA+falpha-2)*fD);
    double cosApam2D = cos((fA+falpha-2)*fD);
    double cosaC12D = cos(falpha*fC12*fD);
    double sinaC12D = sin(falpha*fC12*fD);

    // This is coming straight from the MINOS NueAna package...

    // First we calculate the terms for the alpha expansion (good to all orders in th13)

    // Leading order term 
    double p1 = fcos_sq_th23*fsin_sq_2th13*sinC13D*sinC13D/(fC13*fC13);

    // Terms that appear at order alpha
    double p2Inner =
      fD*cosC13D*(1-fA*fcos_2th13)/fC13 - 
      fA*sinC13D*(fcos_2th13-fA)/(fC13*fC13); 

    double p2 = -2*fsin_sq_th12*fcos_sq_th23*fsin_sq_2th13*sinC13D/(fC13*fC13)*p2Inner*falpha;

    double p3Inner =
      -fsin_dCPproxy*(cosC13D - cos1pAD)*fC13 
      + fcos_dCPproxy*(fC13*sin1pAD - (1-fA*fcos_2th13)*sinC13D);

    double p3 = fsin_2th12*(-fsin_2th23)*fsin_th13*sinC13D/(fA*fC13*fC13)*p3Inner*falpha;

    //  p1 + p2 + p3 is the complete contribution for this expansion
  
    // Now for the expansion in orders of sin(th13) (good to all order falpha) 

    double pa1 = 0.0, pa2 = 0.0;
    if (fabs(falpha)>1E-10) {
      // leading order term
      pa1 = fsin_th23*fsin_th23*fsin_sq_2th12*sinaC12D*sinaC12D/(fC12*fC12);

      // the first order in s13 term
      double t1 = (fcos_2th12 - fA/falpha)/fC12 
	- falpha*fA*fC12*fsin_sq_2th12/(2*(1-falpha)*fC12*fC12);
      double t2 = -fcos_dCPproxy*(sinApam2D-sinaC12D*t1);
      double t3 = -(cosaC12D-cosApam2D)*fsin_dCPproxy;
      double denom = (1-fA-falpha+fA*falpha*fcos_th12*fcos_th12)*fC12;
      double t4 = fsin_2th12*(-fsin_2th23)*(1-falpha)*sinaC12D/denom;

      pa2 = t4*(t3+t2)*fsin_th13;
    }
    // pa1+pa2 is the complete contribution from this expansion

    // Now we need to add the two expansions and subtract off the terms that are
    // in both (falpha^1, s13^1)

    double t1 = sinAD*cosdpD*sinAm1D/(fA*(fA-1));
    double repeated = 2*falpha*fsin_2th12*(-fsin_2th23)*fsin_th13*t1;

    // Calculate the total probability
    double totalP = p1+p2+p3 + (pa1+pa2) - repeated;
    return totalP;
  }

  // --------------------------------------------
  double OscCalculator::P_internal_ee(double E, bool antinu, bool fliptime)
  {
    UpdateBasic();
    UpdateEDep(E,antinu,fliptime);

    double cosC13D = cos(fC13*fD);
    double sinC13D = sin(fC13*fD);
    double sinaC12D = sin(falpha*fC12*fD);

    // This is coming straight from the MINOS NueAna package...

    // First we calculate the terms for the alpha expansion (good to all orders in th13)

    // Leading order term 
    double p1 = 1 - fsin_sq_2th13*sinC13D*sinC13D/(fC13*fC13);

    // Terms that appear at order alpha
    double p2Inner =
      fD*cosC13D*(1-fA*fcos_2th13)/fC13 -
      fA*sinC13D*(fcos_2th13-fA)/(fC13*fC13);

    double p2 = 2*fsin_th12*fsin_th12*fsin_sq_2th13*sinC13D/(fC13*fC13)*p2Inner*falpha;

    //  p1 + p2 is the complete contribution for this expansion

    // Now for the expansion in orders of sin(th13) (good to all order alpha)

    double pa1 = 1.0, pa2 = 0.0;
    if (fabs(falpha)>1E-10) {
      // leading order term
      pa1 = 1 - fsin_sq_2th12*sinaC12D*sinaC12D/(fC12*fC12);
    }
    // pa1 is the complete contribution from this expansion, there is no order s13^1 term

    // Now we need to add the two expansions and subtract off the terms that are
    // in both (falpha^1, s13^1)

    double repeated = 1;

    //  Calculate the total probability
    double totalP = p1+p2 + (pa1+pa2) - repeated;
    return totalP;
  }

  // --------------------------------------------
  double OscCalculator::P_internal_mt(double E, bool antinu, bool fliptime)
  {
    UpdateBasic();
    UpdateEDep(E,antinu,fliptime);

    double cosC13D = cos(fC13*fD);
    double sinC13D = sin(fC13*fD);
    double sin1pAD = sin((fA+1)*fD);
    double cos1pAD = cos((fA+1)*fD);
    double sinAD = sin(fA*fD);
    double sinAm1D = sin((fA-1)*fD);
    double cosAm1D = cos((fA-1)*fD);
    double sinApam2D = sin((fA+falpha-2)*fD);
    double cosApam2D = cos((fA+falpha-2)*fD);
    double cosaC12D = cos(falpha*fC12*fD);
    double sinaC12D = sin(falpha*fC12*fD);
    double sin1pAmCD = sin(0.5*(fA+1-fC13)*fD);
    double sin1pApCD = sin(0.5*(fA+1+fC13)*fD);
    double sinD = sin(fD);
    double sin2D = sin(2*fD);
    double cosaC12pApam2D = cos((falpha*fC12+fA+falpha-2)*fD);

    // This is coming straight from the MINOS NueAna package...

    // First we calculate the terms for the alpha expansion (good to all orders in th13)

    // Leading order term 
    double pmt_0 = 0.5*fsin_sq_2th23;
    pmt_0 *= (1 - (fcos_2th13-fA)/fC13)*sin1pAmCD*sin1pAmCD 
      +  (1 + (fcos_2th13-fA)/fC13)*sin1pApCD*sin1pApCD
      - 0.5*fsin_sq_2th13*sinC13D*sinC13D/(fC13*fC13);

    // Terms that appear at order alpha
    double t0, t1, t2, t3;
    t0 = (fcos_th12*fcos_th12-fsin_th12*fsin_th12*fsin_th13*fsin_th13
          *(1+2*fsin_th13*fsin_th13*fA+fA*fA)/(fC13*fC13))*cosC13D*sin1pAD*2;
    t1 = 2*(fcos_th12*fcos_th12*fcos_th13*fcos_th13-fcos_th12*fcos_th12*fsin_th13*fsin_th13
	    +fsin_th12*fsin_th12*fsin_th13*fsin_th13
	    +(fsin_th12*fsin_th12*fsin_th13*fsin_th13-fcos_th12*fcos_th12)*fA);
    t1 *= sinC13D*cos1pAD/fC13;

    t2 =  fsin_th12*fsin_th12*fsin_sq_2th13*sinC13D/(fC13*fC13*fC13);
    t2 *= fA/fD*sin1pAD+fA/fD*(fcos_2th13-fA)/fC13*sinC13D
      - (1-fA*fcos_2th13)*cosC13D;

    double pmt_1 = -0.5*fsin_sq_2th23*fD*(t0+t1+t2);   

    t0 = cosC13D-cos1pAD;
    t1 = 2*fcos_th13*fcos_th13*fsin_dCPproxy*sinC13D/fC13*t0;
    t2 = -fcos_2th23*fcos_dCPproxy*(1+fA)*t0*t0;

    t3  = fcos_2th23*fcos_dCPproxy*(sin1pAD+(fcos_2th13-fA)/fC13*sinC13D);
    t3 *= (1+2*fsin_th13*fsin_th13*fA + fA*fA)*sinC13D/fC13 - (1+fA)*sin1pAD;

    pmt_1 += (t1+t2+t3)*fsin_th13*fsin_2th12*fsin_2th23/(2*fA*fcos_th13*fcos_th13);
    pmt_1 *= falpha;

    //  pmt_0 + pmt_1 is the complete contribution for this expansion

    // Now for the expansion in orders of sin(th13) (good to all order alpha)

    // Leading order term
    double pmt_a0 =  0.5*fsin_sq_2th23;

    pmt_a0 *= 1 - 0.5*fsin_sq_2th12*sinaC12D*sinaC12D/(fC12*fC12)
      - cosaC12pApam2D
      - (1 - (fcos_2th12 - fA/falpha)/fC12)*sinaC12D*sinApam2D;
            
    double denom = (1-fA-falpha+fA*falpha*fcos_th12*fcos_th12)*fC12;

    t0 = (cosaC12D-cosApam2D)*(cosaC12D-cosApam2D);
    t1 = (fcos_2th12 - fA/falpha)/fC12*sinaC12D+sinApam2D;
    t2 = ((fcos_2th12 - fA/falpha)/fC12+2*(1-falpha)/(falpha*fA*fC12))*sinaC12D + sinApam2D;

    t3 = (falpha*fA*fC12)/2*fcos_2th23*fcos_dCPproxy*(t0 + t1*t2);
    t3 += fsin_dCPproxy*(1-falpha)*(cosaC12D-cosApam2D)*sinaC12D;

    double pmt_a1 = fsin_th13*fsin_2th12*fsin_2th23/denom*t3;

    // pmt_a1+pmt_a2 is the complete contribution from this expansion

    // Now we need to add the two expansions and subtract off the terms that are
    // in both (falpha^1, s13^1)

    t1 = fsin_dCPproxy*sinD*sinAD*sinAm1D/(fA*(fA-1));
    t2 = -1/(fA-1)*fcos_dCPproxy*sinD*(fA*sinD-sinAD*cosAm1D/fA)*fcos_2th23/denom;

    t0 =  2*falpha*fsin_2th12*fsin_2th23*fsin_th13*(t1+t2);

    t1 = fsin_sq_2th23*sinD*sinD 
      - falpha*fsin_sq_2th23*fcos_th12*fcos_th12*fD*sin2D;

    double repeated = t0+t1;

    //  Calculate the total probability
    double totalP = pmt_0 + pmt_1 + pmt_a0 + pmt_a1 - repeated;

    return totalP;
  }
}//namespace
