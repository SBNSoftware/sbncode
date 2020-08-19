#include "HiggsDecay.h"

evgen::ldm::HiggsDecay evgen::ldm::BuildDecay(const HiggsFlux &flux, const simb::MCTruth &mct, float pot) {
  HiggsDecay ret;
  ret.parent_mom = flux.mom;
  ret.parent_start = flux.pos;
  ret.decay = mct.GetParticle(0).Position();
  ret.daughterA_mom = mct.GetParticle(0).Momentum();
  ret.daughterA_pdgid = mct.GetParticle(0).PdgCode();
  ret.daughterB_mom = mct.GetParticle(1).Momentum();
  ret.daughterB_pdgid = mct.GetParticle(1).PdgCode();
  ret.pot = pot;
  return ret;
}
