#include "CAFAna/Core/OscCalcSterileApprox.h"

#include "CAFAna/Core/MathUtil.h"

#include "TMath.h"

#include "Math/SpecFunc.h"

#include <cassert>
#include <iostream>
#include <map>

namespace ana
{
  // --------------------------------------------------------------------------
  double OscCalcSterileApprox::GetSinSq2ThetaEE() const
  {
    // The three angles are coupled via
    //
    // ss2thmm = 4*Um4^2*(1-Um4^2)
    // ss2thee = 4*Ue4^2*(1-Ue4^2)
    // ss2thme = 4*Um4^2*Ue4^2
    //
    // Solve for the final (nue survival angle). Emperically, choosing the
    // negative sign in the expression works well.

    // A couple of limits
    if(fSinSq2ThetaMuMu <= 0) return fSinSq2ThetaMuE / 4;
    if(fSinSq2ThetaMuMu >= 1) return fSinSq2ThetaMuE / 2;

    // General case
    const double Ue4sq = fSinSq2ThetaMuE/(2*fSinSq2ThetaMuMu)*(1-sqrt(1-fSinSq2ThetaMuMu));
    assert(!isinf(Ue4sq) && !isnan(Ue4sq));
    const double sinsq2thetaee = Ue4sq*(1-Ue4sq);

    assert(sinsq2thetaee >= 0 && sinsq2thetaee <= 1);

    // Cross-check that we actually solved the equations correctly
    if(Ue4sq > 0){
      const double Umu4sq = fSinSq2ThetaMuE/(4*Ue4sq);
      assert(fabs(4*Umu4sq*(1-Umu4sq) - fSinSq2ThetaMuMu) < 1e-6);
    }

    return sinsq2thetaee;
  }

  // --------------------------------------------------------------------------
  double Si(double x)
  {
    if(isinf(x)) return TMath::Pi()/2;
    return ROOT::Math::sinint(x);
  }

  // --------------------------------------------------------------------------
  double AvgSinSq(double k, double a, double b)
  {
    // This function shows up high in profiles, and is often called with the same masses/baselines/energies
    static std::map<std::tuple<double, double, double>, double> cache;
    const std::tuple<double, double, double> key = {k, a, b};

    // energies times masses times baselines plus some slack
    if(cache.size() > 100*100*10) cache.clear();

    auto it = cache.find(key);
    if(it != cache.end()) return it->second;

    double ret = 0;
    if(a == b){
      ret = util::sqr(sin(k/a));
    }
    else{
      // https://www.wolframalpha.com/input/?i=integral+sin%5E2(k%2Fx)+from+a+to+b
      ret = k*(Si(2*k/a)-Si(2*k/b));
      assert(!isnan(ret));
      if(a) ret -= a * util::sqr(sin(k/a));
      assert(!isnan(ret));
      if(b) ret += b * util::sqr(sin(k/b));
      assert(!isnan(ret));

      ret /= (b-a);
    }

    cache[key] = ret;
    return ret;
  }

  // --------------------------------------------------------------------------
  double OscCalcSterileApprox::P(int from, int to, double E)
  {
    return P_range(from, to, E, E);
  }

  // --------------------------------------------------------------------------
  double OscCalcSterileApprox::PFromDelta(int from, int to, double Delta) const
  {
    if(abs(from) == 14 && abs(to) == 14){
      return 1-fSinSq2ThetaMuMu*Delta;
    }
    else if(abs(from) == 12 && abs(to) == 12){
      return 1-GetSinSq2ThetaEE()*Delta;
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

    //Option to return the active fraction
    else if (abs(from) == 14 && to == 0) {
      return (1-fSinSq2ThetaMuMu*Delta) + (fSinSq2ThetaMuE*Delta);
    }
    else if (abs(from) == 12 && to == 0) {
      return (1-GetSinSq2ThetaEE()*Delta) + (fSinSq2ThetaMuE*Delta);
    }

    std::cout << "OscCalcSterileApprox: P(" << from << ", " << to << ") not implemented" << std::endl;
    abort();
  }

  // --------------------------------------------------------------------------
  double OscCalcSterileApprox::P_range(int from, int to, double Elo, double Ehi)
  {
    if(Ehi <= 0) return 0;

    if(fDmsq == 0) return (from == to || to == 0) ? 1 : 0;

    Elo = std::max(0., Elo);

    assert(fL != 0);
    const double Delta = AvgSinSq(1.267*fDmsq*fL, Elo, Ehi);
    assert(!isnan(Delta));

    return PFromDelta(from, to, Delta);
  }

  // --------------------------------------------------------------------------
  double OscCalcSterileApprox::P_LoverE(int from, int to,
                                        double LElo, double LEhi)
  {
    if(fDmsq == 0) return (from == to || to == 0) ? 1 : 0;

    LElo = std::max(0., LElo);

    LElo *= 1.267*fDmsq;
    LEhi *= 1.267*fDmsq;

    if(LElo == LEhi) return PFromDelta(from, to, util::sqr(sin(LElo)));

    // Average value of sin^2 over the range
    const double Delta = (.25*(sin(2*LElo) - sin(2*LEhi)) + .5*(LEhi - LElo))/(LEhi-LElo);
    assert(!isnan(Delta));

    return PFromDelta(from, to, Delta);
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
  TMD5* OscCalcSterileApprox::GetParamsHash() const
  {
    TMD5* ret = new TMD5;
    const std::string txt = "SterileApprox";
    ret->Update((unsigned char*)txt.c_str(), txt.size());
    const int kNumParams = 4;
    double buf[kNumParams] = {fDmsq, fSinSq2ThetaMuMu, fSinSq2ThetaMuE, fL};
    ret->Update((unsigned char*)buf, sizeof(double)*kNumParams);
    ret->Final();
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
  const OscCalcSterileApprox* DowncastToSterileApprox(const osc::IOscCalc* calc, bool allowFail)
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
  OscCalcSterileApprox* DowncastToSterileApprox(osc::IOscCalc* calc,
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
