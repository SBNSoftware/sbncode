#ifndef BASICTOOL_GEODIRECTEDLINE_CXX
#define BASICTOOL_GEODIRECTEDLINE_CXX

#include "GeoDirectedLine.h"

namespace geoalgo {

  DirectedLine::DirectedLine() : Line()
  {}

  DirectedLine::DirectedLine(const double x, const double y, const double z,
			     const double dirx, const double diry, const double dirz)
    : Line( x, y, z, x+dirx, y+diry, z+dirz)
  { check_and_raise(_pt1,_pt2); }
    
  DirectedLine::DirectedLine(const Point_t& pt, const Vector_t& dir)
    : Line( pt, pt+dir)
  { check_and_raise(_pt1,_pt2); }

  DirectedLine::DirectedLine(const HalfLine& l)
    : Line( l.Start(), l.Start()+l.Dir() )
  { check_and_raise(_pt1,_pt2); }

  Vector_t DirectedLine::Dir() const
  { return _pt2 - _pt1; }
  
}
#endif


