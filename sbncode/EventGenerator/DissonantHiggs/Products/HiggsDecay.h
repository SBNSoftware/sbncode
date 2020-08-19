#ifndef _HiggsDecay_HH_
#define _HiggsDecay_HH_

#include "TLorentzVector.h"
#include "HiggsFlux.h"
#include "nusimdata/SimulationBase/MCTruth.h"

namespace evgen {
namespace ldm {
class HiggsDecay {
public:
  TLorentzVector parent_mom;
  TLorentzVector parent_start;
  TLorentzVector decay;
  TLorentzVector daughterA_mom;
  float daughterA_pdgid;
  TLorentzVector daughterB_mom;
  float daughterB_pdgid;
  float pot;
};

HiggsDecay BuildDecay(const HiggsFlux &flux, const simb::MCTruth &mct, float pot);

} // end namespace ldm

} // end namespace evgen

#endif
