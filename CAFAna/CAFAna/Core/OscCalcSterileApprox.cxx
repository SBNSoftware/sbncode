#include "CAFAna/Core/OscCalcSterileApprox.h"

#include "Utilities/func/MathUtil.h"

#include "TMath.h"

#include "Math/SpecFunc.h"

#include <cassert>
#include <iostream>

namespace ana
{
  // --------------------------------------------------------------------------
  double Si(double x)
  {
    if(isinf(x)) return TMath::Pi()/2;
    return ROOT::Math::sinint(x);
  }

  // --------------------------------------------------------------------------
  double AvgSinSq(double k, double a, double b)
  {
    if(a == b) return util::sqr(sin(k));

    // https://www.wolframalpha.com/input/?i=integral+sin%5E2(k%2Fx)+from+a+to+b
    double ret = k*(Si(2*k/a)-Si(2*k/b));
    assert(!isnan(ret));
    if(a) ret -= a * util::sqr(sin(k/a));
    assert(!isnan(ret));
    if(b) ret += b * util::sqr(sin(k/b));
    assert(!isnan(ret));

    ret /= (b-a);

    return ret;
  }

  // --------------------------------------------------------------------------
  double OscCalcSterileApprox::P(int from, int to, double E)
  {
    return P(from, to, E, E);
  }

  // --------------------------------------------------------------------------
  double OscCalcSterileApprox::P(int from, int to, double Elo, double Ehi)
  {
    if(Ehi <= 0) return 0;

    if(fDmsq == 0) return (from == to) ? 1 : 0;

    Elo = std::max(0., Elo);

    assert(fL != 0);
    const double Delta = AvgSinSq(1.267*fDmsq*fL, Elo, Ehi);
    assert(!isnan(Delta));

    // The three angles are coupled via
    //
    // ss2thmm = 4*Um4^2*(1-Um4^2)
    // ss2thee = 4*Ue4^2*(1-Ue4^2)
    // ss2thme = 4*Um4^2*Ue4^2
    //
    // Solve for the final (nue survival angle), and take the smaller angle
    // solution
    // TODO TODO TODO
    const double sinsq2thetaee = 0/*fSinSq2ThetaMuMu > 0 ? fSinSq2ThetaMuE * (1 - sqrt(1-4*fSinSq2ThetaMuMu))/(2*fSinSq2ThetaMuMu) : 0*/;

    if(abs(from) == 14 && abs(to) == 14){
      return 1-fSinSq2ThetaMuMu*Delta;
    }
    else if(abs(from) == 12 && abs(to) == 12){
      return 1-sinsq2thetaee*Delta;
    }
    else if(abs(from) == 14 && abs(to) == 12){
      return fSinSq2ThetaMuE*Delta;
    }
    else if(abs(from) == 12 && abs(to) == 14){
      // TODO - this seems reasonable, is it right?
      return fSinSq2ThetaMuE*Delta;
    }
    else if(abs(to) == 16){ // no tau appearance
      return 0;
    }

    std::cout << "OscCalculatorSterileApprox: P(" << from << ", " << to << ") not implemented" << std::endl;
    abort();
  }

  // --------------------------------------------------------------------------
  OscCalcSterileApprox* OscCalcSterileApprox::Copy() const
  {
    OscCalcSterileApprox* ret = new OscCalcSterileApprox;

    ret->fDmsq = fDmsq;
    ret->fSinSq2ThetaMuMu = fSinSq2ThetaMuMu;
    ret->fSinSq2ThetaMuE = fSinSq2ThetaMuE;
    ret->fL = fL;

    return ret;
  }

  //---------------------------------------------------------------------------
  OscCalcSterileApproxAdjustable* DefaultSterileApproxCalc()
  {
    auto ret = new OscCalcSterileApproxAdjustable;
    ret->calc.SetDmsq(0);
    ret->calc.SetSinSq2ThetaMuMu(0);
    ret->calc.SetSinSq2ThetaMuE(0);
    ret->calc.SetL(0); // make clear this is uninitialized

    return ret;
  }

  //---------------------------------------------------------------------------
  const OscCalcSterileApprox* DowncastToSterileApprox(const osc::IOscCalculator* calc, bool allowFail)
  {
    const OscCalcSterileApprox* ret1
      = dynamic_cast<const OscCalcSterileApprox*>(calc);
    if(ret1) return ret1;

    const OscCalcSterileApproxAdjustable* ret2
      = dynamic_cast<const OscCalcSterileApproxAdjustable*>(calc);
    if(ret2) return &ret2->calc;

    if(allowFail) return 0;

    std::cout << "Calculator was not of type OscCalcSterileApprox." << std::endl;
    abort();
  }

  //---------------------------------------------------------------------------
  OscCalcSterileApprox* DowncastToSterileApprox(osc::IOscCalculator* calc,
                                                bool allowFail)
  {
    OscCalcSterileApprox* ret1
      = dynamic_cast<OscCalcSterileApprox*>(calc);
    if(ret1) return ret1;

    OscCalcSterileApproxAdjustable* ret2
      = dynamic_cast<OscCalcSterileApproxAdjustable*>(calc);
    if(ret2) return &ret2->calc;

    if(allowFail) return 0;

    std::cout << "Calculator was not of type OscCalcSterileApprox." << std::endl;
    abort();
  }

}
