#include "DissonantHiggs.h"

evgen::ldm::DissonantHiggs evgen::ldm::BuildHiggs(const HiggsFlux &flux, const HiggsDecay &decay, std::array<TVector3, 2> inout, double flux_weight, double ray_weight, double decay_weight, double pot) {
  DissonantHiggs ret;
  ret.kaon_dmom = flux.kmom;
  ret.kaon_dmom_beamcoord = flux.kmom_beamcoord;
  ret.kaon_pdg = flux.kaon_pdg;
  ret.kaon_dpos_beamcoord = flux.pos_beamcoord;

  ret.higgs_mom = flux.mom;
  ret.higgs_start = flux.pos;
  ret.higgs_mom_beamcoord = flux.mom_beamcoord;

  ret.decay_pos = decay.pos;
  ret.daughterA_mom = decay.daughterA_mom;
  ret.daughterA_pdg = decay.daughterA_pdg;
  ret.daughterB_mom = decay.daughterB_mom;
  ret.daughterB_pdg = decay.daughterB_pdg;

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
