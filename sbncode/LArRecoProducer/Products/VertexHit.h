#ifndef sbncode_VertexHit_HH
#define sbncode_VertexHit_HH

#include <vector>
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

namespace sbn {
  class VertexHit {
  public:
    geo::WireID wire;    
    float charge;
    float proj_dist_to_vertex;
    std::vector<int> nearbyPFPIDs;
    std::vector<float> PCproj;
  };
} // end namespace sbn

#endif
