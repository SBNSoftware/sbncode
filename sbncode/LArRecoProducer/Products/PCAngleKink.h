#ifndef sbncode_PCAngleKink_HH
#define sbncode_PCAngleKink_HH

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

namespace sbn {
  class PCAngleKink {
  public:
    geo::WireID maxWire;
    float fwhm_distance;
    float est_angle;
    float max_angle;

    float fit_chi2;
    float fit_angle;
    float fit_pitch;

    std::array<float, 2> position_max;
    std::array<float, 2> position_lo;
    std::array<float, 2> position_hi;

    std::array<float, 2> vec_lo_at_max;
    std::array<float, 2> vec_hi_at_max;
    std::array<float, 2> vec_lo_at_halfmax_lo;
    std::array<float, 2> vec_hi_at_halfmax_lo;
    std::array<float, 2> vec_lo_at_halfmax_hi;
    std::array<float, 2> vec_hi_at_halfmax_hi;
  };
} // end namespace sbn

#endif
