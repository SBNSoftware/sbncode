#include "CAFAna/Systs/SystComponentScale.h"

#include "TObjString.h"
#include "TDirectory.h"

#include <cmath>
#include <iostream>

namespace ana
{
  //----------------------------------------------------------------------
  SystComponentScale::~SystComponentScale()
  {
  }

  //----------------------------------------------------------------------
  void SystComponentScale::Shift(double sigma,
                                 caf::SRSliceProxy* slc,
                                 double& weight) const
  {
    if(!fCut(slc)) return;

    if(fType == kExponential){
      weight *= pow(1+fOneSigma, sigma);
    }
    else{
      weight *= 1+sigma*fOneSigma;
      weight = std::max(0., weight);
    }
  }

  //----------------------------------------------------------------------
  std::unique_ptr<SystComponentScale> SystComponentScale::LoadFrom(TDirectory* dir)
  {
    TObjString* ptag = (TObjString*)dir->Get("type");
    assert(ptag);

    const TString tag = ptag->GetString();

    std::cerr << "Unknown SystComponentScale type '" << tag << "'" << std::endl;
    abort();
  }
}
