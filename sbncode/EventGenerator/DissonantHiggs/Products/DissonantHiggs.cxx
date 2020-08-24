#include "DissonantHiggs.h"

evgen::ldm::DissonantHiggs evgen::ldm::BuildHiggs(const HiggsFlux &flux, const KaonParent &kaon, const HiggsDecay &decay, const simb::MCTruth &mct, std::array<TVector3, 2> inout, double flux_weight, double ray_weight, double decay_weight, double pot) {
  DissonantHiggs ret;
  ret.kaon_dmom = kaon.mom;
  ret.kaon_dpos = kaon.pos;
  ret.kaon_pdg = kaon.kaon_pdg;
  ret.kaon_wgt = kaon.weight;
  ret.higgs_mom = flux.mom;
  ret.higgs_start = flux.pos;
  ret.higgs_mom_beamcoord = flux.mom_beamcoord;

  if (mct.NParticles() > 0) {
    ret.decay_pos = mct.GetParticle(0).Position();
    ret.daughterA_mom = mct.GetParticle(0).Momentum();
    ret.daughterA_pdgid = mct.GetParticle(0).PdgCode();
    ret.daughterB_mom = mct.GetParticle(1).Momentum();
    ret.daughterB_pdgid = mct.GetParticle(1).PdgCode();
  }

  ret.pot = pot;
  ret.flux_weight = flux_weight;
  ret.ray_weight = ray_weight;
  ret.decay_weight = decay_weight;
  ret.higgs_enter = inout[0];
  ret.higgs_exit = inout[1];

  ret.mass = flux.mass;
  ret.mixing = flux.mixing;

  ret.decay_width = decay.decay_width;
  ret.mean_lifetime = decay.mean_lifetime;
  ret.mean_distance = decay.mean_distance;

  return ret;
}
