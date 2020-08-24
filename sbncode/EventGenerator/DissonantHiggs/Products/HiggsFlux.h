#ifndef _HiggsFlux_HH_
#define _HiggsFlux_HH_

#include "TLorentzVector.h"

namespace evgen {
namespace ldm {
class HiggsFlux {
public:
  TLorentzVector pos;
  TLorentzVector kmom;
  TLorentzVector mom;
  TLorentzVector mom_beamcoord;
  float mixing;
  float mass;

  double decay_width;
  double mean_lifetime;
  double mean_distance;
};
}
}

#endif
