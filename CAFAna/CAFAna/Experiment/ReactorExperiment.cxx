#include "CAFAna/Experiment/ReactorExperiment.h"

#include "CAFAna/Vars/FitVars.h"

#include "Utilities/func/MathUtil.h"

#include "TDirectory.h"
#include "TH1.h"
#include "TObjString.h"

#include <cassert>

namespace ana
{
  //----------------------------------------------------------------------
  double ReactorExperiment::ChiSq(osc::IOscCalculatorAdjustable* osc,
                                  const SystShifts& /*syst*/) const
  {
    return util::sqr((kFitSinSq2Theta13.GetValue(osc)-fBestFit)/fSigma);
  }

  //----------------------------------------------------------------------
  void ReactorExperiment::SaveTo(TDirectory* dir) const
  {
    TDirectory* tmp = dir;

    dir->cd();
    TObjString("ReactorExperiment").Write("type");

    TH1D params("", "", 2, 0, 2);
    params.SetBinContent(1, fBestFit);
    params.SetBinContent(2, fSigma);
    params.Write("params");

    tmp->cd();
  }

  //----------------------------------------------------------------------
  std::unique_ptr<ReactorExperiment> ReactorExperiment::LoadFrom(TDirectory* dir)
  {
    TObjString* tag = (TObjString*)dir->Get("type");
    assert(tag);
    assert(tag->GetString() == "ReactorExperiment");

    TH1* params = (TH1*)dir->Get("params");
    assert(params);

    const double bestFit = params->GetBinContent(1);
    const double sigma   = params->GetBinContent(2);

    return std::unique_ptr<ReactorExperiment>(new ReactorExperiment(bestFit,
                                                                    sigma));
  }

  //----------------------------------------------------------------------
  const ReactorExperiment* DayaBayConstraint2014()
  {
    // Daya Bay, Neutrino 2014:
    // https://indico.fnal.gov/getFile.py/access?contribId=256&sessionId=15&resId=0&materialId=slides&confId=8022
    return new ReactorExperiment(.084, .005);
  }

  //----------------------------------------------------------------------
  const ReactorExperiment* WorldReactorConstraint2015()
  {
    // Weighted average of
    // Daya Bay: arXiv:1505.03456
    // RENO: arXiv:1410.7987
    // Double Chooz: 1406.7763
    return new ReactorExperiment(.086, .005);
  }

  //----------------------------------------------------------------------
  const ReactorExperiment* WorldReactorConstraint2016()
  {
    // PDG website as of Jun 2016
    return new ReactorExperiment(.085, .005);
  }

  //----------------------------------------------------------------------
  const ReactorExperiment* WorldReactorConstraint2017()
  {
    // http://pdg.lbl.gov/2017/tables/rpp2017-sum-leptons.pdf
    // ssth13=0.021+/-0.0011 -> ss2th13=0.082+/-0.004
    return new ReactorExperiment(.082, .004);
  }
}
