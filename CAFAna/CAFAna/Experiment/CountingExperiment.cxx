#include "CAFAna/Experiment/CountingExperiment.h"

#include "CAFAna/Prediction/IPrediction.h"

#include "CAFAna/Core/LoadFromFile.h"
#include "CAFAna/Core/Utilities.h"

#include "OscLib/func/IOscCalculator.h"

#include "TDirectory.h"
#include "TObjString.h"
#include "TH1.h"

namespace ana
{
  //----------------------------------------------------------------------
  CountingExperiment::CountingExperiment(const IPrediction* p,
                                         const Spectrum& d,
                                         const Spectrum& cosmic)
    : fMC(p), fData(d),
      fCosmic(cosmic.ToTH1(d.Livetime(), kLivetime))
  {
  }

  //----------------------------------------------------------------------
  CountingExperiment::~CountingExperiment()
  {
    delete fCosmic;
  }

  //----------------------------------------------------------------------
  double CountingExperiment::ChiSq(osc::IOscCalculatorAdjustable* osc,
                                   const SystShifts& syst) const
  {
    double exp = fMC->PredictSyst(osc, syst).Integral(fData.POT());
    if (fCosmic) exp += fCosmic->Integral(0,-1);
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
    fData.SaveTo(dir->mkdir("data"));
    
    if(fCosmic) fCosmic->Write("cosmic");

    tmp->cd();
  }

  //----------------------------------------------------------------------
  std::unique_ptr<CountingExperiment> CountingExperiment::LoadFrom(TDirectory* dir)
  {
    TObjString* ptag = (TObjString*)dir->Get("type");
    assert(ptag);
    assert(ptag->GetString() == "CountingExperiment");

    assert(dir->GetDirectory("mc"));
    assert(dir->GetDirectory("data"));
    

    const IPrediction* mc = ana::LoadFrom<IPrediction>(dir->GetDirectory("mc")).release();
    const std::unique_ptr<Spectrum> data = Spectrum::LoadFrom(dir->GetDirectory("data"));

    TH1* cosmic = 0;
    if(dir->Get("cosmic")) cosmic = (TH1*)dir->Get("cosmic");

    auto ret = std::make_unique<CountingExperiment>(mc, *data);
    if(cosmic) ret->fCosmic = cosmic;
    return ret;
  }
}
