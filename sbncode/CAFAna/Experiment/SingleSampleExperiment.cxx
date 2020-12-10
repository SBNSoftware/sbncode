#include "CAFAna/Experiment/SingleSampleExperiment.h"

#include "CAFAna/Core/LoadFromFile.h"
//#include "CAFAna/Core/StanUtils.h"
#include "CAFAna/Core/Utilities.h"

#include "OscLib/IOscCalc.h"

#include "TDirectory.h"
#include "TObjString.h"
#include "TH2.h"

namespace ana
{
  //----------------------------------------------------------------------
  SingleSampleExperiment::SingleSampleExperiment(const IPrediction* pred,
                                                 const Spectrum& data)
    : fMC(pred), fData(data)
  {
  }

  //----------------------------------------------------------------------
  SingleSampleExperiment::~SingleSampleExperiment()
  {
  }

  //----------------------------------------------------------------------
  double SingleSampleExperiment::ChiSq(osc::IOscCalcAdjustable* calc,
                                       const SystShifts& syst) const
  {
    Eigen::ArrayXd apred = fMC->PredictSyst(calc, syst).GetEigen(fData.POT());
    Eigen::ArrayXd adata = fData.GetEigen(fData.POT());

    ApplyMask(apred, adata);

    // full namespace qualification to avoid degeneracy with method inherited
    // from IExperiment
    return ana::LogLikelihood(apred, adata);
  }

  //----------------------------------------------------------------------
  void SingleSampleExperiment::ApplyMask(Eigen::ArrayXd& a,
                                         Eigen::ArrayXd& b) const
  {
    if(fMaskA.size() == 0) return;

    assert(a.size() == fMaskA.size());
    assert(b.size() == fMaskA.size());

    // Arrays mean we get bin-by-bin operations
    a *= fMaskA;
    b *= fMaskA;
  }

  //----------------------------------------------------------------------
  /*
  stan::math::var SingleSampleExperiment::LogLikelihood(osc::IOscCalcAdjustableStan *osc,
                                                        const SystShifts &syst) const
  {
    const Spectrum pred = fMC->PredictSyst(osc, syst);

    const Eigen::ArrayXd data = fData.GetEigen(fData.POT());

    // It's possible to have a non-stan prediction. e.g. from a NoOsc
    // prediction with no systs.
    if(pred.HasStan()){
      // fully-qualified so that we get the one in StanUtils.h
      //
      // LogLikelihood(), confusingly, returns chi2=-2*LL
      return ana::LogLikelihood(pred.GetEigenStan(fData.POT()), data) / -2.;
    }
    else{
      return ana::LogLikelihood(pred.GetEigen(fData.POT()), data) / -2.;
    }
  }
  */

  //----------------------------------------------------------------------
  void SingleSampleExperiment::SaveTo(TDirectory* dir) const
  {
    TDirectory* tmp = dir;

    dir->cd();
    TObjString("SingleSampleExperiment").Write("type");

    fMC->SaveTo(dir->mkdir("mc"));
    fData.SaveTo(dir, "data");

    tmp->cd();
  }

  //----------------------------------------------------------------------
  std::unique_ptr<SingleSampleExperiment> SingleSampleExperiment::LoadFrom(TDirectory* dir)
  {
    TObjString* ptag = (TObjString*)dir->Get("type");
    assert(ptag);
    assert(ptag->GetString() == "SingleSampleExperiment");

    assert(dir->GetDirectory("mc"));
    

    const IPrediction* mc = ana::LoadFrom<IPrediction>(dir->GetDirectory("mc")).release();
    const std::unique_ptr<Spectrum> data = Spectrum::LoadFrom(dir, "data");

    auto ret = std::make_unique<SingleSampleExperiment>(mc, *data);
    return ret;
  }

  //----------------------------------------------------------------------
  void SingleSampleExperiment::SetMaskHist(double xmin, double xmax, double ymin, double ymax)
  {
    fMaskA = GetMaskArray(fData, xmin, xmax, ymin, ymax);
  }
}
