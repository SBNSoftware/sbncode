#ifndef BASICTOOL_GEOAABOX_CXX
#define BASICTOOL_GEOAABOX_CXX

#include "GeoAABox.h"

namespace geoalgo {

  AABox::AABox() 
    : _min(3)
    , _max(3)
  {}
    
  AABox::AABox(const double x_min, const double y_min, const double z_min,
	       const double x_max, const double y_max, const double z_max)
    : _min ( x_min, y_min, z_min )
    , _max ( x_max, y_max, z_max )
  {}
    
  AABox::AABox(const Point_t& min, const Vector_t& max)
    : _min ( min )
    , _max ( max   )
  { 
    if(min.size()!=3 || max.size()!=3)
      throw GeoAlgoException("AABox ctor accepts only 3D Point!");
  }
  
  const Point_t& AABox::Min() const { return _min; }
  const Point_t& AABox::Max() const { return _max; }

  void AABox::Min(const double x, const double y, const double z)
  { _min[0] = x; _min[1] = y; _min[2] = z; }
  void AABox::Max(const double x, const double y, const double z)
  { _max[0] = x; _max[1] = y; _max[2] = z; }

  bool AABox::Contain(const Point_t &pt) const {
    return !( (pt[0] < _min[0] || _max[0] < pt[0]) || // point is outside X boundaries OR
	      (pt[1] < _min[1] || _max[1] < pt[1]) || // point is outside Y boundaries OR
	      (pt[2] < _min[2] || _max[2] < pt[2])    // point is outside Z boundaries
	      );
  }

}
#endif


