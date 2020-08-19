#ifndef _HiggsFlux_HH_
#define _HiggsFlux_HH_

#include "TLorentzVector.h"

namespace evgen {
namespace ldm {
class HiggsFlux {
public:
  TLorentzVector pos;
  TLorentzVector mom;
};
}
}

#endif
