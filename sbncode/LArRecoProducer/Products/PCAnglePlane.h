#ifndef sbncode_PCAnglePlane_HH
#define sbncode_PCAnglePlane_HH

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

namespace sbn {
  class PCAngle {
  public:
    geo::WireID wire;
    std::array<float, 2> lo_vector;
    std::array<float, 2> hi_vector;
    float angle;
    float hit_time_cm;
    float hit_wire_cm;
    float dist_to_lo;
    float dist_to_hi;
    int hitID;
    int hitIDLo;
    int hitIDHi;
    bool complete;
  };

  class PCAnglePlane {
  public:
    std::vector<std::vector<PCAngle>> angles;
    std::vector<std::vector<int>> branchHierarchy;
    std::vector<int> branchIDs;
    std::vector<int> generations;
    geo::PlaneID plane;
    unsigned nBranches;
  };

} // end namespace sbn

#endif
