#pragma once

#include "CAFAna/Core/ISyst.h"

#include "TString.h"

class TH1;

namespace ana {

class SBNFluxHadronSyst : public ISyst
{
public:
  virtual ~SBNFluxHadronSyst();

  void Shift(double sigma, caf::SRSliceProxy* slc, double& weight) const override;

protected:
  friend const SBNFluxHadronSyst* GetSBNFluxHadronSyst(unsigned int);

  SBNFluxHadronSyst(int i)
      : ISyst(TString::Format("flux_%i", i).Data(),
              TString::Format("FluxHadron #%i", i).Data()),
        fIdx(i), fScale() {}

  int fIdx;

  mutable TH1* fScale[3];
};

const SBNFluxHadronSyst* GetSBNFluxHadronSyst(unsigned int i);

// Because vector<T*> won't automatically convert to vector<U*> even when U
// inherits from V.
struct SBNFluxHadronSystVector : public std::vector<const SBNFluxHadronSyst*> {
  operator std::vector<const ISyst*>() {
    return std::vector<const ISyst*>(begin(), end());
  }
};

/// \param N total number of systematics
SBNFluxHadronSystVector GetSBNFluxHadronSysts(unsigned int N);

} // namespace ana
