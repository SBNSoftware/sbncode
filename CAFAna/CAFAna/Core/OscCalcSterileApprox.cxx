#include "CAFAna/Core/OscCalcSterileApprox.h"

#include "Utilities/func/MathUtil.h"

#include "TMath.h"

#include "Math/SpecFunc.h"

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
    // https://www.wolframalpha.com/input/?i=integral+sin%5E2(k%2Fx)+from+a+to+b
    double ret = k*(Si(2*k/a)-Si(2*k/b));
    if(a) ret -= a * util::sqr(sin(k/a));
    if(b) ret += b * util::sqr(sin(k/b));

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
    Elo = std::max(0., Elo);

    const double arg = 1.267*fDmsq*fL;
    const double Delta = (Elo != Ehi) ? AvgSinSq(arg, Elo, Ehi) : util::sqr(sin(arg));

    // The three angles are coupled via
    //
    // ss2thmm = 4*Um4^2*(1-Um4^2)
    // ss2thee = 4*Ue4^2*(1-Ue4^2)
    // ss2thme = 4*Um4^2*Ue4^2
    //
    // Solve for the final (nue survival angle), and take the smaller angle
    // solution
    const double sinsq2thetaee = fSinSq2ThetaMuMu > 0 ? fSinSq2ThetaMuE * (1 - sqrt(1-4*fSinSq2ThetaMuMu))/(2*fSinSq2ThetaMuMu) : 0;

    if(abs(from) == 14 && abs(to) == 14){
      return 1-fSinSq2ThetaMuMu*Delta;
    }
    else if(abs(from) == 12 && abs(to) == 12){
      return 1-sinsq2thetaee*Delta;
    }
    else if(abs(from) == 14 && abs(to) == 12){
      return fSinSq2ThetaMuE*Delta;
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
}
