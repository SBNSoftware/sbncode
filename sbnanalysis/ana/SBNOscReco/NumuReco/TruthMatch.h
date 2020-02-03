#ifndef _sbncode_TruthMatch_hh
#define _sbncode_TruthMatch_hh

#include <vector>
#include "../Data/RecoEvent.h"
#include "core/Event.hh"

// Do truth match with the NumuReco objects defined in Data/

namespace numu {
  void ApplySliceTruthMatch(RecoEvent &event, const std::vector<event::Interaction> &truth);
}



#endif
