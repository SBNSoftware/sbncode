#ifndef _DissonantHiggs_HH_
#define _DissonantHiggs_HH_

#include "TLorentzVector.h"
#include "HiggsDecay.h"
#include "HiggsFlux.h"
#include <array>

namespace evgen {
namespace ldm {
class DissonantHiggs {
public:
  TLorentzVector kaon_dmom;
  TLorentzVector kaon_dmom_beamcoord;
  TLorentzVector kaon_dpos_beamcoord;
  int kaon_pdg;
  TLorentzVector higgs_mom_beamcoord;
  TLorentzVector higgs_mom;
  TLorentzVector higgs_start;
  TVector3 higgs_enter;
  TVector3 higgs_exit;
  TLorentzVector decay_pos;
  TLorentzVector daughterA_mom;
  int daughterA_pdg;
  TLorentzVector daughterB_mom;
  int daughterB_pdg;
  double pot;
  double flux_weight;
  double ray_weight;
  double decay_weight;

  double mass;
  double mixing;

  double decay_width;
  double mean_lifetime;
  double mean_distance;
};

DissonantHiggs BuildHiggs(const HiggsFlux &flux, const HiggsDecay &decay, std::array<TVector3, 2> inout, double flux_weight, double ray_weight, double decay_weight, double pot);

} // end namespace ldm

} // end namespace evgen

#endif
