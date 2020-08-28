#ifndef _HiggsDecay_HH_
#define _HiggsDecay_HH_

#include "TLorentzVector.h"

namespace evgen {
namespace ldm {
class HiggsDecay {
public:
  TLorentzVector pos;
  TLorentzVector daughterA_mom;
  int daughterA_pdg;
  TLorentzVector daughterB_mom;
  int daughterB_pdg;
  double daughter_mass;
  double decay_width;
  double mean_lifetime;
  double mean_distance;
};

} // end namespace ldm

} // end namespace evgen

#endif
