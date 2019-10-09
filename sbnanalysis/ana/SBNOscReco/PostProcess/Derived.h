#ifndef _sbncode_Derived_hh
#define _sbncode_Derived_hh

#include <vector>
#include "../Data/RecoEvent.h"

// helper functions to calculate derived quantities from
// output classes of the reconstruction processing

namespace numu {
  float dist2Match(const numu::RecoInteraction &vertex, const std::vector<numu::RecoInteraction> &candidates);
  float trackMatchCompletion(unsigned truth_index, const numu::RecoEvent &event);
}



#endif
