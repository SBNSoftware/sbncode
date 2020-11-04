#ifndef BASICTOOL_GEOLINESEGMENT_CXX
#define BASICTOOL_GEOLINESEGMENT_CXX

#include "GeoLineSegment.h"

namespace geoalgo {

  LineSegment::LineSegment()
    : _start(3)
    , _end(3)
    , _dir(3)
  {DirReset();}
    
  LineSegment::LineSegment(const double start_x, const double start_y, const double start_z,
			   const double end_x,   const double end_y,   const double end_z   )
    : _start ( start_x, start_y, start_z )
    , _end   ( end_x,   end_y,   end_z   )
    , _dir   (3)
  {DirReset();}

  LineSegment::LineSegment(const Point_t& start, const Point_t& end)
    : _start ( start )
    , _end   ( end   )
    , _dir   (3)
  { 
    if(start.size()!=3 || end.size()!=3)
      throw GeoAlgoException("LineSegment ctor accepts only 3D Point!");
    DirReset();
  }

  const Point_t& LineSegment::Start() const { return _start; }

  const Point_t& LineSegment::End() const { return _end; }

  const Vector_t LineSegment::Dir() const { return _dir; }

  void LineSegment::Start(const double x, const double y, const double z)
  { _start[0] = x; _start[1] = y; _start[2] = z; 
    DirReset();
  }

  void LineSegment::End(const double x, const double y, const double z)
  { _end[0] = x; _end[1] = y; _end[2] = z; 
    DirReset();
  }
  
  void LineSegment::DirReset() { _dir = _end - _start; }

}

#endif

