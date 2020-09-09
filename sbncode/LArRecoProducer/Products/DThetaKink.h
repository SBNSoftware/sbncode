#ifndef sbncode_DThetaKink_HH
#define sbncode_DThetaKink_HH

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

namespace sbn {
  class DThetaKink {
    public:
      geo::WireID wire;
      float hit_time_cm;
      float hit_wire_cm;
      float pca;
      float vec;

      bool complete;
      unsigned branch;
      unsigned generation;

      int hitID;
      int hitIDLo;
      int hitIDHi;
  };
} // end namespace sbn

#endif
