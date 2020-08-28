#ifndef _HiggsFlux_HH_
#define _HiggsFlux_HH_

#include "TLorentzVector.h"

namespace evgen {
namespace ldm {
class HiggsFlux {
public:
  TLorentzVector pos_beamcoord;
  TLorentzVector pos;
  TLorentzVector kmom_beamcoord;
  TLorentzVector kmom;
  TLorentzVector mom;
  TLorentzVector mom_beamcoord;
  float mixing;
  float mass;
  int kaon_pdg;

};
}
}

#endif
