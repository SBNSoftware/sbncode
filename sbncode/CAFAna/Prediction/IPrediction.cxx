#include "CAFAna/Prediction/IPrediction.h"

#include "CAFAna/Core/LoadFromFile.h"

#include "OscLib/func/IOscCalculator.h"

#include <cassert>
#include <iostream>

#include "TDirectory.h"
#include "TObjString.h"

// To implement LoadFrom()
#include "CAFAna/Prediction/PredictionNoExtrap.h"
#include "CAFAna/Prediction/PredictionInterp.h"
#include "CAFAna/Prediction/PredictionNoOsc.h"
#include "CAFAna/Prediction/PredictionScaleComp.h"

namespace ana
{
  //----------------------------------------------------------------------
  // Definition to satisfy declaration in Core/LoadFromFile.h
  template<> std::unique_ptr<IPrediction> LoadFrom<IPrediction>(TDirectory* dir)
  {
    TObjString* ptag = (TObjString*)dir->Get("type");
    assert(ptag);

    const TString tag = ptag->GetString();

    if(tag == "PredictionNoExtrap") return PredictionNoExtrap::LoadFrom(dir);

    // Backwards compatibility
    if(tag == "PredictionInterp" ||
       tag == "PredictionInterp2") return PredictionInterp::LoadFrom(dir);

    if(tag == "PredictionNoOsc") return PredictionNoOsc::LoadFrom(dir);

    if(tag == "PredictionScaleComp") return PredictionScaleComp::LoadFrom(dir);

    std::cerr << "Unknown Prediction type '" << tag << "'" << std::endl;
    abort();
  }

  //----------------------------------------------------------------------
  Spectrum IPrediction::PredictUnoscillated() const
  {
    // Default implementation
    osc::NoOscillations noosc;
    return Predict(&noosc);
  }

  //----------------------------------------------------------------------
  Spectrum IPrediction::PredictSyst(osc::IOscCalculator* calc,
                                    const SystShifts& syst) const
  {
    assert(syst.IsNominal() && "This Prediction doesn't support PredictSyst(). Did you just mean Predict()?");

    // Default implementation: no treatment of systematics
    return Predict(calc);
  }

  //----------------------------------------------------------------------
  Spectrum IPrediction::PredictComponentSyst(osc::IOscCalculator* calc,
                                             const SystShifts& syst,
                                             Flavors::Flavors_t flav,
                                             Current::Current_t curr,
                                             Sign::Sign_t sign) const
  {
    assert(syst.IsNominal() && "This Prediction doesn't support PredictSyst(). Did you just mean Predict()?");

    // Default implementation: no treatment of systematics
    return PredictComponent(calc, flav, curr, sign);
  }

  //----------------------------------------------------------------------
  void IPrediction::SaveTo(TDirectory* dir) const
  {
    assert(0 && "Not implemented");
  }
}
