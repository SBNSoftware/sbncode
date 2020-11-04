#pragma once

#include "OscLib/IOscCalc.h"

namespace ana
{
  class OscCalcSterileApprox: public osc::IOscCalc
  {
  public:
    // if flavAfter == 0, give the active fraction
    virtual double P(int from, int to, double E) override;
    double P(int from, int to, double Elo, double Ehi);

    double P_LoverE(int from, int to, double LElo, double LEhi);

    virtual OscCalcSterileApprox* Copy() const override;

    void SetDmsq(double d) {fDmsq = d;}
    double GetDmsq() const {return fDmsq;}

    void SetSinSq2ThetaMuMu(double t) {fSinSq2ThetaMuMu = t;}
    double GetSinSq2ThetaMuMu() const {return fSinSq2ThetaMuMu;}

    void SetSinSq2ThetaMuE(double t) {fSinSq2ThetaMuE = t;}
    double GetSinSq2ThetaMuE() const {return fSinSq2ThetaMuE;}

    double GetSinSq2ThetaEE() const; ///< calculated from the others

    // TODO - potentially remove L in the brave new L/E future
    void SetL(double L) {fL = L;}
    double GetL() const {return fL;}

    TMD5* GetParamsHash() const override;
  protected:
    double PFromDelta(int from, int to, double Delta) const;

    double fDmsq;
    double fSinSq2ThetaMuMu;
    double fSinSq2ThetaMuE;
    double fL;
  };

  class OscCalcSterileApproxAdjustable: public osc::IOscCalcAdjustable
  {
  public:
    OscCalcSterileApprox calc;

    virtual double P(int from, int to, double E) override
    {
      return calc.P(from, to, E);
    }

    virtual OscCalcSterileApproxAdjustable* Copy() const override
    {
      auto ret = new OscCalcSterileApproxAdjustable;
      auto c = calc.Copy();
      ret->calc = *c;
      delete c;
      return ret;
    }

    virtual void SetL     (double L     ) override {calc.SetL(L);}
    virtual void SetRho   (double rho   ) override {}
    virtual void SetDmsq21(const double& dmsq21) override {}
    virtual void SetDmsq32(const double& dmsq32) override {}
    virtual void SetTh12  (const double& th12  ) override {}
    virtual void SetTh13  (const double& th13  ) override {}
    virtual void SetTh23  (const double& th23  ) override {}
    virtual void SetdCP   (const double& dCP   ) override {}

    TMD5* GetParamsHash() const override {return calc.GetParamsHash();}
  };

  OscCalcSterileApproxAdjustable* DefaultSterileApproxCalc();

  const OscCalcSterileApprox* DowncastToSterileApprox(const osc::IOscCalc* calc, bool allowFail = false);
  OscCalcSterileApprox* DowncastToSterileApprox(osc::IOscCalc* calc, bool allowFail = false);
}
