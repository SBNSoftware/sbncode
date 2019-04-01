#include "OscLib/func/IOscCalculator.h"

namespace ana
{
  class OscCalcSterileApprox: public osc::IOscCalculator
  {
  public:
    virtual double P(int from, int to, double E) override;
    double P(int from, int to, double Elo, double Ehi);

    virtual OscCalcSterileApprox* Copy() const override;

    void SetDmsq(double d) {fDmsq = d;}
    double GetDmsq() const {return fDmsq;}

    void SetSinSq2ThetaMuMu(double t) {fSinSq2ThetaMuMu = t;}
    double GetSinSq2ThetaMuMu() const {return fSinSq2ThetaMuMu;}

    void SetSinSq2ThetaMuE(double t) {fSinSq2ThetaMuE = t;}
    double GetSinSq2ThetaMuE() const {return fSinSq2ThetaMuE;}

    void SetL(double L) {fL = L;}
    double GetL() const {return fL;}
  protected:
    double fDmsq;
    double fSinSq2ThetaMuMu;
    double fSinSq2ThetaMuE;
    double fL;
  };

  class OscCalcSterileApproxAdjustable: osc::IOscCalculatorAdjustable
  {
  public:
    OscCalcSterileApprox calc;
  };
}
