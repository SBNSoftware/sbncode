#ifndef sbncode_MergedTrack_HH
#define sbncode_MergedTrack_HH

#include "TVector3.h"

namespace sbn {
  class MergedTrackInfo {
  public:
    // std::array<bool, 3> trunk_wire_direction_is_ascending;
    // std::array<int, 3> branch_wire_start;
    TVector3 vertex;
    TVector3 direction;
    int trunk;
    int branch;
    float branch_overlap;
    float branch_start;
    float trunk_start;
  };
} // end namespace sbn

#endif
