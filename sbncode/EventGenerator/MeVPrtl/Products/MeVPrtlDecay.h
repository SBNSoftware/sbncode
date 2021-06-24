#ifndef _MeVPrtlDecay_HH_
#define _MeVPrtlDecay_HH_

#include "TLorentzVector.h"

namespace evgen {
namespace ldm {
class MeVPrtlDecay {
public:
  TLorentzVector pos;
  std::vector<TLorentzVector> daughter_mom;
  std::vector<int> daughter_pdg;

  double decay_width;
  double mean_lifetime;
  double mean_distance;
};

} // end namespace ldm

} // end namespace evgen

#endif
