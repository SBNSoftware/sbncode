#include "CAFAna/Analysis/Calcs.h"

#include "OscLib/func/OscCalculatorPMNSOpt.h"
#include "OscLib/func/OscCalculatorSterile.h"

#include <cmath>
#include <iostream>

namespace ana
{
  //----------------------------------------------------------------------
  void ResetOscCalcToDefault(osc::IOscCalculatorAdjustable* calc)
  {
    // Default to ICARUS
    calc->SetL(0.6);

    // TODO: get good reference from Alex and Joao
    calc->SetRho(2.84); // g/cm^3

    // PDG: http://pdg.lbl.gov/2014/tables/rpp2014-sum-leptons.pdf
    calc->SetDmsq21(7.53e-5);
    calc->SetTh12(asin(sqrt(.846))/2);

    // Picking anything other than maximal mixing opens a can of worms
    calc->SetTh23(M_PI/4);

    // NH from PDG2014 + 2015 update
    //http://pdg.lbl.gov/2015/tables/rpp2015-sum-leptons.pdf 
    calc->SetDmsq32(2.44e-3); // NB: this is normal hierarchy

    // Reactor average from PDG2014 + 2015 update
    calc->SetTh13(asin(sqrt(.085))/2);

    // Going to have to plot for nue analysis anyway
    calc->SetdCP(0);
  }

  //----------------------------------------------------------------------
  osc::IOscCalculatorAdjustable* DefaultOscCalc()
  {
    osc::IOscCalculatorAdjustable* ret = new osc::OscCalculatorPMNSOpt;
    ResetOscCalcToDefault(ret);
    return ret;
  }

  //----------------------------------------------------------------------
  void ResetOscCalcToDefaultIH(osc::IOscCalculatorAdjustable* calc)
  {
    //Share most defaults 
    ResetOscCalcToDefault(calc);
    // IH from PDG2014 + 2015 update 
    // http://pdg.lbl.gov/2015/tables/rpp2015-sum-leptons.pdf
    calc->SetDmsq32(-2.49e-3); 

  }

  //----------------------------------------------------------------------
  osc::IOscCalculatorAdjustable* DefaultOscCalcIH()
  {
    osc::IOscCalculatorAdjustable* ret = new osc::OscCalculatorPMNSOpt;
    ResetOscCalcToDefaultIH(ret);
    return ret;
  }

  //----------------------------------------------------------------------
  void ResetSterileCalcToDefault(osc::OscCalculatorSterile* calc)
  {
    osc::OscCalculatorPMNSOpt* tmp = new osc::OscCalculatorPMNSOpt();
    ResetOscCalcToDefault(tmp);

    calc->SetL(tmp->GetL());
    calc->SetRho(tmp->GetRho());

    calc->SetDm(2, tmp->GetDmsq21());
    calc->SetDm(3, tmp->GetDmsq21() + tmp->GetDmsq32());

    calc->SetAngle(1, 2, tmp->GetTh12());
    calc->SetAngle(1, 3, tmp->GetTh13());
    calc->SetAngle(2, 3, tmp->GetTh23());

    calc->SetDelta(1, 3, tmp->GetdCP());

    delete tmp;
  }

  //----------------------------------------------------------------------
  osc::OscCalculatorSterile* DefaultSterileCalc(int nflavors)
  {
    osc::OscCalculatorSterile* ret = new osc::OscCalculatorSterile;

    if(nflavors < 3) {
      std::cout << "The default calculator requires at least 3 flavors." << std::endl;
      std::cout << "Using 3 flavors." << std::endl;
      ret->SetNFlavors(3);
    }
    else {
      ret->SetNFlavors(nflavors);
    }

    ResetSterileCalcToDefault(ret);
    return ret;
  }
}
