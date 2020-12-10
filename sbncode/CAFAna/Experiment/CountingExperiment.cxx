#include "CAFAna/Experiment/CountingExperiment.h"

#include "CAFAna/Prediction/IPrediction.h"

#include "CAFAna/Core/LoadFromFile.h"
#include "CAFAna/Core/Utilities.h"

#include "OscLib/IOscCalc.h"

#include "TDirectory.h"
#include "TObjString.h"

namespace ana
{
  //----------------------------------------------------------------------
  CountingExperiment::~CountingExperiment()
  {
  }

  //----------------------------------------------------------------------
  double CountingExperiment::ChiSq(osc::IOscCalcAdjustable* osc,
                                   const SystShifts& syst) const
  {
    double exp = fMC->PredictSyst(osc, syst).Integral(fData.POT());
    double obs = fData.Integral(fData.POT());

    return LogLikelihood(exp, obs);
  }

  //----------------------------------------------------------------------
  void CountingExperiment::SaveTo(TDirectory* dir) const
  {
    TDirectory* tmp = dir;

    dir->cd();
    TObjString("CountingExperiment").Write("type");

    fMC->SaveTo(dir->mkdir("mc"));
    fData.SaveTo(dir, "data");

    tmp->cd();
  }

  //----------------------------------------------------------------------
  std::unique_ptr<CountingExperiment> CountingExperiment::LoadFrom(TDirectory* dir)
  {
    TObjString* ptag = (TObjString*)dir->Get("type");
    assert(ptag);
    assert(ptag->GetString() == "CountingExperiment");

    assert(dir->GetDirectory("mc"));
    

    const IPrediction* mc = ana::LoadFrom<IPrediction>(dir->GetDirectory("mc")).release();
    const std::unique_ptr<Spectrum> data = Spectrum::LoadFrom(dir, "data");

    auto ret = std::make_unique<CountingExperiment>(mc, *data);
    return ret;
  }
}
