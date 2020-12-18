#include "MeVPrtlTruth.h"

evgen::ldm::MeVPrtlTruth evgen::ldm::BuildMeVPrtlTruth(const MeVPrtlFlux &flux, const MeVPrtlDecay &decay, std::array<TVector3, 2> inout, double flux_weight, double ray_weight, double decay_weight, double pot) {
  MeVPrtlTruth ret;
  ret.kaon_dmom = flux.kmom;
  ret.kaon_dmom_beamcoord = flux.kmom_beamcoord;
  ret.kaon_pdg = flux.kaon_pdg;
  ret.kaon_dpos_beamcoord = flux.pos_beamcoord;

  ret.mevprtl_mom = flux.mom;
  ret.mevprtl_start = flux.pos;
  ret.mevprtl_mom_beamcoord = flux.mom_beamcoord;

  ret.decay_pos = decay.pos;
  ret.daughterA_mom = decay.daughterA_mom;
  ret.daughterA_pdg = decay.daughterA_pdg;
  ret.daughterB_mom = decay.daughterB_mom;
  ret.daughterB_pdg = decay.daughterB_pdg;

  ret.pot = pot;
  ret.flux_weight = flux_weight;
  ret.ray_weight = ray_weight;
  ret.decay_weight = decay_weight;
  ret.mevprtl_enter = inout[0];
  ret.mevprtl_exit = inout[1];

  ret.mass = flux.mass;
  ret.C1 = flux.C1;
  ret.C2 = flux.C2;
  ret.C3 = flux.C3;
  ret.C4 = flux.C4;
  ret.C5 = flux.C5;

  ret.decay_width = decay.decay_width;
  ret.mean_lifetime = decay.mean_lifetime;
  ret.mean_distance = decay.mean_distance;

  ret.gen = (evgen::ldm::Generator)flux.generator;

  return ret;
}
