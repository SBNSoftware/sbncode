#include "MeVPrtlTruth.h"

evgen::ldm::MeVPrtlTruth::MeVPrtlTruth(const MeVPrtlFlux &flux, const MeVPrtlDecay &decay, std::array<TVector3, 2> inout, double flux_weight, double ray_weight, double decay_weight, double pot) {
  kaon_dmom = flux.kmom;
  kaon_dmom_beamcoord = flux.kmom_beamcoord;
  kaon_pdg = flux.kaon_pdg;
  kaon_dpos_beamcoord = flux.pos_beamcoord;

  mevprtl_mom = flux.mom;
  mevprtl_start = flux.pos;
  mevprtl_mom_beamcoord = flux.mom_beamcoord;

  decay_pos = decay.pos;

  daughter_mom = decay.daughter_mom;
  daughter_pdg = decay.daughter_pdg;

  pot = pot;
  flux_weight = flux_weight;
  ray_weight = ray_weight;
  decay_weight = decay_weight;
  mevprtl_enter = inout[0];
  mevprtl_exit = inout[1];

  mass = flux.mass;
  C1 = flux.C1;
  C2 = flux.C2;
  C3 = flux.C3;
  C4 = flux.C4;
  C5 = flux.C5;

  decay_width = decay.decay_width;
  mean_lifetime = decay.mean_lifetime;
  mean_distance = decay.mean_distance;

  gen = (evgen::ldm::Generator)flux.generator;
}
