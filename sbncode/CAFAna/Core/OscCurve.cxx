#include "CAFAna/Core/OscCurve.h"

#include "CAFAna/Core/Binning.h"

#include "OscLib/IOscCalc.h"

#include "TH1.h"

#include <iostream>
#include <map>

#include "CAFAna/Core/OscCalcSterileApprox.h"

namespace
{
  template<class T> inline bool IsNoOscillations(const osc::_IOscCalc<T>* c)
  {
    return dynamic_cast<const osc::_NoOscillations<T>*>(c) != 0;
  }
}

namespace ana
{
  //----------------------------------------------------------------------
  /// Helper for constructors
  template<class T> Eigen::Array<T, Eigen::Dynamic, 1>
  ToEigen(osc::_IOscCalc<T>* calc, int from, int to)
  {
    const unsigned int N = kTrueLOverEBins.NBins();

    // Have to allow for underflow and overflow
    Eigen::Array<T, Eigen::Dynamic, 1> ret(N+2);
    ret[0] = (from == to || to == 0) ? 1 : 0; // underflow
    ret[N+1] = 0; // overflow


    // We have extra knowledge that calculators of this type have special
    // modes allowing calculation in L/E and an intrinsic energy smearing.
    OscCalcSterileApprox* approx = DowncastToSterileApprox(calc, true);

    if(approx){
      for(unsigned int i = 1; i <= N; ++i){
        const double LElo = kTrueLOverEBins.Edges()[i-1];
        const double LEhi = kTrueLOverEBins.Edges()[i];
        ret[i] = approx->P_LoverE(from, to, LElo, LEhi);
      }
    }
    else{
      if(!IsNoOscillations(calc)){
        std::cout << "Trying to use a calculator which is not OscCalcSterileApprox with an L/E axis. Will have to code up additional hacks for this to work" << std::endl;
        abort();
      }

      for(unsigned int i = 1; i <= N; ++i) ret[i] = calc->P(from, to, 1); // energy irrelevant for NoOsc calc
    }

    return ret;
  }

  //----------------------------------------------------------------------
  OscCurve::OscCurve(osc::IOscCalc* calc, int from, int to)
    : Ratio(Hist::Adopt(ToEigen(calc, from, to)),
            std::vector<Binning>(1, kTrueLOverEBins),
            std::vector<std::string>(1, "True L / E (km / GeV)")),
      fFrom(from), fTo(to)
  {
  }

  //----------------------------------------------------------------------
  OscCurve::OscCurve(osc::IOscCalcStan* calc, int from, int to)
    : Ratio(Hist::AdoptStan(ToEigen(calc, from, to)),
            std::vector<Binning>(1, kTrueLOverEBins),
            std::vector<std::string>(1, "True L / E (km / GeV)")),
      fFrom(from), fTo(to)
  {
  }

  //----------------------------------------------------------------------
  OscCurve::~OscCurve()
  {
  }

  //----------------------------------------------------------------------
  TH1D* OscCurve::ToTH1(bool title) const
  {
    // Could have a file temporarily open
    DontAddDirectory guard;

    TH1D* ret = Ratio::ToTH1();
    ret->GetYaxis()->SetTitle("Probability");

    if(title){
      // Don't do this work unless it's explicitly requested
      std::map<int, std::string> nus;
      nus[12] = nus[-12] = "e";
      nus[14] = nus[-14] = "#mu";
      nus[16] = nus[-16] = "#tau";
      nus[0] = "active";
      const std::string nu = (fFrom > 0) ? "#nu" : "#bar#nu";

      ret->SetTitle((nu+"_{"+nus[fFrom]+"}#rightarrow"+nu+"_{"+nus[fTo]+"}").c_str());
    }

    return ret;
  }
}
