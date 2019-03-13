#include "CAFAna/Systs/SystComponentScale.h"

// #include "CAFAna/Systs/DUNEXSecSysts.h"

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
                                  Restorer& /*restore*/,
                                  caf::StandardRecord* sr,
                                  double& weight) const
  {
    if(!fCut(sr)) return;

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

    // if(tag == "DUNEXSecSyst") return DUNEXSecSyst::LoadFrom(dir);

    std::cerr << "Unknown SystComponentScale type '" << tag << "'" << std::endl;
    abort();
  }
}
