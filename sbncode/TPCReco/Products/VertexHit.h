#ifndef sbncode_VertexHit_HH
#define sbncode_VertexHit_HH

#include <vector>
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "TVector3.h"

namespace sbn {
  class VertexHit {
  public:
    geo::WireID wire;    
    float charge;
    float proj_dist_to_vertex;
    int spID;
    TVector3 spXYZ;
    float pitch;
    float dqdx;
    float dedx;
    std::vector<int> nearbyPFPIDs;
    std::vector<float> nearbyPFPDists;
    std::vector<float> nearbyPFP3DDists;
    std::vector<int> nearbyPFP3DIDs;
  };
} // end namespace sbn

#endif
