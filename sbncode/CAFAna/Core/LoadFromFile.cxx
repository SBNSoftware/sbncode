#include "CAFAna/Core/LoadFromFile.h"

#include "CAFAna/Core/OscCalcSterileApprox.h"

#include "OscLib/OscCalc.h"
#include "OscLib/OscCalcDumb.h"
#include "OscLib/OscCalcGeneral.h"
#include "OscLib/OscCalcPMNS.h"
#include "OscLib/OscCalcPMNSOpt.h"
#include "OscLib/OscCalcSterile.h"
#include "OscLib/OscCalcPMNS_NSI.h"

#include "TObjString.h"
#include "TH1.h"
#include "TVectorD.h"
#include <vector>

namespace ana
{
  //----------------------------------------------------------------------
  template<> std::unique_ptr<osc::IOscCalc>
  LoadFrom<osc::IOscCalc>(TDirectory* dir)
  {
    TObjString* ptag = (TObjString*)dir->Get("type");
    assert(ptag);
    const TString tag = ptag->GetString();

    if(tag == "NoOscillations") return std::unique_ptr<osc::IOscCalc>(new osc::NoOscillations);
    if(tag == "OscCalculatorDumb") return std::unique_ptr<osc::IOscCalc>(new osc::OscCalcDumb);

    TVectorD* params = (TVectorD*)dir->Get("params");
    assert(params);

    osc::IOscCalcAdjustable* ret = 0;

    if(tag == "OscCalculator") ret = new osc::OscCalc;
    if(tag == "OscCalculatorGeneral") ret = new osc::OscCalcGeneral;
    if(tag == "OscCalculatorPMNS") ret = new osc::OscCalcPMNS;
    if(tag == "OscCalculatorPMNSOpt") ret = new osc::OscCalcPMNSOpt;
    if(tag == "OscCalculatorSterile") ret = new osc::OscCalcSterile;
    if(tag == "OscCalculatorPMNS_NSI") ret = new osc::OscCalcPMNS_NSI;
    if(tag == "OscCalcSterileApprox"){
      assert(params->GetNrows() == 4);
      auto ret = new OscCalcSterileApproxAdjustable;
      ret->calc.SetDmsq((*params)[0]);
      ret->calc.SetSinSq2ThetaMuMu((*params)[1]);
      ret->calc.SetSinSq2ThetaMuE((*params)[2]);
      ret->calc.SetL((*params)[3]);
      return std::unique_ptr<osc::IOscCalcAdjustable>(ret);
    }

    if(!ret){
      std::cout << "LoadFrom not implemented for " << tag << std::endl;
      abort();
    }

    //special case how OscCalcSterile is initialized
    if(tag == "OscCalculatorSterile") {
      std::vector<double> state;
      for (int i = 0; i < params->GetNrows(); ++i) {
	state.push_back( (*params)[i] );
      }
      dynamic_cast<osc::OscCalcSterile*>(ret)->SetState(state);
      return std::unique_ptr<osc::IOscCalcAdjustable>(ret);
    }
    //special case how OscCalcPMNS_NSI is initialized
    if(tag == "OscCalculatorPMNS_NSI") {
      std::vector<double> state;
      for (int i = 0; i < params->GetNrows(); ++i) {
	state.push_back( (*params)[i] );
      }
      dynamic_cast<osc::OscCalcPMNS_NSI*>(ret)->SetState(state);
      return std::unique_ptr<osc::IOscCalcAdjustable>(ret);
    }

    assert(params->GetNrows() == 8);

    ret->SetL     ((*params)[0]);
    ret->SetRho   ((*params)[1]);
    ret->SetDmsq21((*params)[2]);
    ret->SetDmsq32((*params)[3]);
    ret->SetTh12  ((*params)[4]);
    ret->SetTh13  ((*params)[5]);
    ret->SetTh23  ((*params)[6]);
    ret->SetdCP   ((*params)[7]);

    return std::unique_ptr<osc::IOscCalcAdjustable>(ret);
  }

  //----------------------------------------------------------------------
  template<> void SaveTo(const osc::IOscCalc& x, TDirectory* dir)
  {
    TDirectory* tmp = gDirectory;

    dir->cd();

    if(dynamic_cast<const osc::NoOscillations*>(&x)){
      TObjString("NoOscillations").Write("type");
      tmp->cd();
      return;
    }

    if(dynamic_cast<const osc::OscCalcDumb*>(&x)){
      TObjString("OscCalculatorDumb").Write("type");
      tmp->cd();
      return;
    }

    const OscCalcSterileApproxAdjustable* approx = dynamic_cast<const OscCalcSterileApproxAdjustable*>(&x);
    if(approx){
      TObjString("OscCalcSterileApprox").Write("type");
      TVectorD params(4);
      params[0] = approx->calc.GetDmsq();
      params[1] = approx->calc.GetSinSq2ThetaMuMu();
      params[2] = approx->calc.GetSinSq2ThetaMuE();
      params[3] = approx->calc.GetL();
      params.Write("params");
      tmp->cd();
      return;
    }

    const osc::IOscCalcAdjustable* y = dynamic_cast<const osc::IOscCalcAdjustable*>(&x);
    if(!y){
      std::cout << "Unknown calculator in SaveTo " << typeid(x).name() << std::endl;
      abort();
    }
    
    const osc::OscCalcSterile* tmpSterile = dynamic_cast<const osc::OscCalcSterile*>(&x);
    if(tmpSterile) {
      TObjString("OscCalculatorSterile").Write("type");      
      std::vector<double> state = tmpSterile->GetState();
      TVectorD params(state.size());
      for (unsigned int i = 0; i < state.size(); ++i) params[i] = state[i];
      params.Write("params");
      tmp->cd();
      return;
    }

    //for the implementation of OscCalcPMNS_NSI
    const osc::OscCalcPMNS_NSI* tmpNSI = dynamic_cast<const osc::OscCalcPMNS_NSI*>(&x);
    if(tmpNSI) {
      TObjString("OscCalcPMNS_NSI").Write("type");      
      std::vector<double> state = tmpNSI->GetState();
      TVectorD params(state.size());
      for (unsigned int i = 0; i < state.size(); ++i) params[i] = state[i];
      params.Write("params");
      tmp->cd();
      return;
    }

    /* */if(dynamic_cast<const osc::OscCalc*>(&x)) TObjString("OscCalculatorPMNS").Write("type");
    else if(dynamic_cast<const osc::OscCalcGeneral*>(&x)) TObjString("OscCalculatorGeneral").Write("type");
    else if(dynamic_cast<const osc::OscCalcPMNS*>(&x)) TObjString("OscCalculatorPMNS").Write("type");
    else if(dynamic_cast<const osc::OscCalcPMNSOpt*>(&x)) TObjString("OscCalculatorPMNSOpt").Write("type");
    else{
      std::cout << "Unimplemented calculator in SaveTo " << typeid(x).name() << std::endl;
      abort();
    }
    
    TVectorD params(8);

    params[0] = y->GetL();
    params[1] = y->GetRho();
    params[2] = y->GetDmsq21();
    params[3] = y->GetDmsq32();
    params[4] = y->GetTh12();
    params[5] = y->GetTh13();
    params[6] = y->GetTh23();
    params[7] = y->GetdCP();

    params.Write("params");

    tmp->cd();
  }
}
