#ifndef _DissonantHiggs_HH_
#define _DissonantHiggs_HH_

#include "TLorentzVector.h"
#include "HiggsDecay.h"
#include "HiggsFlux.h"
#include "KaonParent.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include <array>

namespace evgen {
namespace ldm {
class DissonantHiggs {
public:
  TLorentzVector kaon_dmom;
  TLorentzVector kaon_dpos;
  int kaon_pdg;
  double kaon_wgt;
  TLorentzVector higgs_mom_beamcoord;
  TLorentzVector higgs_mom;
  TLorentzVector higgs_start;
  TVector3 higgs_enter;
  TVector3 higgs_exit;
  TLorentzVector decay_pos;
  TLorentzVector daughterA_mom;
  float daughterA_pdgid;
  TLorentzVector daughterB_mom;
  float daughterB_pdgid;
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

DissonantHiggs BuildHiggs(const HiggsFlux &flux, const KaonParent &kaon, const HiggsDecay &decay, const simb::MCTruth &mct, std::array<TVector3, 2> inout, double flux_weight, double ray_weight, double decay_weight, double pot);

} // end namespace ldm

} // end namespace evgen

#endif
