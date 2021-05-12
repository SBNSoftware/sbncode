#pragma once

#include "CAFAna/Core/ISyst.h"
#include "StandardRecord/Proxy/SRProxy.h"
#include "CAFAna/Core/Utilities.h"

#include "TFile.h"
#include "TH1.h"

namespace ana
{
  class MECSyst: public ISyst
  {
  public:
  MECSyst() : ISyst("mec", "MEC Syst") {}
    void Shift(double sigma,
	       caf::SRProxy* sr, double& weight) const override
    {
      if(sigma < -1) sigma = -1;
      if(sr->truth[0].neutrino.genie_intcode == 10) weight *= 1 + sigma;
    }
  };

  const MECSyst& GetMECSyst()
  {
    static const MECSyst mSyst;
    return mSyst;
  }

  class POTSyst: public ISyst
  {
  public:
  POTSyst() : ISyst("potsyst", "2% POT Syst") {}
    void Shift(double sigma,
	       caf::SRProxy* sr, double& weight) const override
    {
      weight *= (1.0 + 0.02*sigma);
    }
  };

  const POTSyst& GetPOTSyst()
  {
    static const POTSyst pSyst;
    return pSyst;
  }
}
