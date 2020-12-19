#ifndef _MeVPrtlFlux_HH_
#define _MeVPrtlFlux_HH_

#include "TLorentzVector.h"

namespace evgen {
namespace ldm {
class MeVPrtlFlux {
public:
  TLorentzVector pos_beamcoord;
  TLorentzVector pos;
  TLorentzVector kmom_beamcoord;
  TLorentzVector kmom;
  TLorentzVector mom;
  TLorentzVector mom_beamcoord;
  float C1;
  float C2;
  float C3;
  float C4;
  float C5;
  float mass;
  int kaon_pdg;
  int secondary_pdg;
  int generator;

};
}
}

#endif
