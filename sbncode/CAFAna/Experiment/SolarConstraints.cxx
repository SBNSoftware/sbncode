#include "CAFAna/Experiment/SolarConstraints.h"

#include "OscLib/func/IOscCalculator.h"
#include "Utilities/func/MathUtil.h"

#include "TDirectory.h"
#include "TH1.h"
#include "TObjString.h"

#include <cassert>

namespace ana
{
  //----------------------------------------------------------------------
  SolarConstraints::SolarConstraints()
  {
    // These value are from the 2014 PDG
    // http://pdg.lbl.gov/2014/tables/rpp2014-sum-leptons.pdf
    std::cerr << "WARNING: Using 2014 Solar constraints."
	      << "Are you sure you don't want kSolarConstraintsPDG2017 ?" 
	      << std::endl;

    fCentralDmsq = 7.53e-5;
    fErrorDmsq   = 0.18e-5;

    fCentralAngle = 0.846;
    fErrorAngle   = 0.021;
  }

  //----------------------------------------------------------------------
  SolarConstraints::SolarConstraints(const double dmsq,  const double errorDmsq,
				     const double ss2th, const double errorSs2th)
    :fCentralDmsq (dmsq),  fErrorDmsq (errorDmsq),
     fCentralAngle(ss2th), fErrorAngle(errorSs2th)
  { }

  //----------------------------------------------------------------------
  double SolarConstraints::ChiSq(osc::IOscCalculatorAdjustable* osc,
                                 const SystShifts& /*syst*/) const
  {
    double ret = 0;

    ret += util::sqr((osc->GetDmsq21() - fCentralDmsq)/fErrorDmsq);

    const double ss2th12 = util::sqr(sin(2*osc->GetTh12()));

    ret += util::sqr((ss2th12 - fCentralAngle)/fErrorAngle);

    return ret;
  }

  //----------------------------------------------------------------------
  void SolarConstraints::SaveTo(TDirectory* dir) const
  {
    TDirectory* tmp = dir;

    dir->cd();
    TObjString("SolarConstraints").Write("type");

    TH1D params("", "", 4, 0, 4);
    params.SetBinContent(1, fCentralDmsq);
    params.SetBinContent(2, fErrorDmsq);
    params.SetBinContent(3, fCentralAngle);
    params.SetBinContent(4, fErrorAngle);
    params.Write("params");

    tmp->cd();
  }

  //----------------------------------------------------------------------
  std::unique_ptr<SolarConstraints> SolarConstraints::LoadFrom(TDirectory* dir)
  {
    TObjString* tag = (TObjString*)dir->Get("type");
    assert(tag);
    assert(tag->GetString() == "SolarConstraints");

    std::unique_ptr<SolarConstraints> ret(new SolarConstraints);

    TH1* params = (TH1*)dir->Get("params");
    assert(params);

    ret->fCentralDmsq  = params->GetBinContent(1);
    ret->fErrorDmsq    = params->GetBinContent(2);
    ret->fCentralAngle = params->GetBinContent(3);
    ret->fErrorAngle   = params->GetBinContent(4);

    return ret;
  }
}
