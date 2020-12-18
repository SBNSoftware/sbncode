#ifndef _MeVPrtlTruth_HH_
#define _MeVPrtlTruth_HH_

#include "TLorentzVector.h"
#include "MeVPrtlDecay.h"
#include "MeVPrtlFlux.h"
#include <array>

namespace evgen {
namespace ldm {

enum Generator {
  kUnknown=-1,
  kDissonantHiggs=0,
  kHNL=1
};

class MeVPrtlTruth {
public:
  TLorentzVector kaon_dmom;
  TLorentzVector kaon_dmom_beamcoord;
  TLorentzVector kaon_dpos_beamcoord;
  int kaon_pdg;
  TLorentzVector mevprtl_mom_beamcoord;
  TLorentzVector mevprtl_mom;
  TLorentzVector mevprtl_start;
  TVector3 mevprtl_enter;
  TVector3 mevprtl_exit;
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
  double C1;
  double C2;
  double C3;
  double C4;
  double C5;

  double decay_width;
  double mean_lifetime;
  double mean_distance;

  Generator gen;
};

MeVPrtlTruth BuildMeVPrtlTruth(const MeVPrtlFlux &flux, const MeVPrtlDecay &decay, std::array<TVector3, 2> inout, double flux_weight, double ray_weight, double decay_weight, double pot);

} // end namespace ldm

} // end namespace evgen

#endif
