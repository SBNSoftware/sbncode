#ifndef _KaonParent_HH_
#define _KaonParent_HH_

#include "TLorentzVector.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "TRotation.h"

namespace evgen {
namespace ldm {
class KaonParent {
public:
  TLorentzVector pos;
  TLorentzVector mom;
  int kaon_pdg;
  int pion_pdg;
  double weight;
};

bool MakeKaonParent(const simb::MCFlux &flux, KaonParent &ret);

}
}

#endif
