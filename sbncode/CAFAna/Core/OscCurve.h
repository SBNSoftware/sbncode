#pragma once

#include <map>
#include <string>

class TH1;
class TH1D;

namespace osc
{
  template<class T> class _IOscCalc;
  typedef _IOscCalc<double> IOscCalc;
}

namespace ana
{
  /// Transition probability for any one channel as a function of energy
  class OscCurve
  {
  public:
    OscCurve(osc::IOscCalc* calc, int from, int to, bool LoverE);
    OscCurve(TH1* h);
    virtual ~OscCurve();

    OscCurve(const OscCurve& rhs);
    OscCurve& operator=(const OscCurve& rhs);

    TH1D* ToTH1(bool title = false) const;
  protected:
    int fFrom, fTo;
    TH1D* fHist;
  };
}
