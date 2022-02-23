#ifndef CAF_FOLDEMDAUGHTERS_H
#define CAF_FOLDEMDAUGHTERS_H

#include "sbnanaobj/StandardRecord/SRTrueParticle.h"

#include <vector>

namespace caf {
  std::vector<caf::SRTrueParticle> FoldEMShowerDaughters(const std::vector<caf::SRTrueParticle> &ps);
}

#endif
