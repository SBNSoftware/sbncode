#ifndef BASICTOOL_GEOLINE_CXX
#define BASICTOOL_GEOLINE_CXX

#include "GeoLine.h"

namespace geoalgo {

  Line::Line() 
    : _pt1(3)
    , _pt2(3)
  {}

  Line::Line(const double x1, const double y1, const double z1,
	     const double x2, const double y2, const double z2)
    : _pt1 (x1, y1, z1)
    , _pt2 (x2, y2, z2)
  { check_and_raise(_pt1,_pt2); }

  Line::Line(const Point_t& pt1, const Point_t& pt2)
    : _pt1 ( pt1 )
    , _pt2 ( pt2 )
  { check_and_raise(pt1,pt2); }
  
  const Point_t& Line::Pt1() const { return _pt1; }
  const Point_t& Line::Pt2() const { return _pt2; }
    
  void Line::Pt1(const double x, const double y, const double z)
  { 
    _pt1[0] = x; 
    _pt1[1] = y; 
    _pt1[2] = z; 
    check_and_raise(_pt1,_pt2);
  }

  void Line::Pt2(const double x, const double y, const double z)
  { 
    _pt2[0] = x; 
    _pt2[1] = y; 
    _pt2[2] = z; 
    check_and_raise(_pt1,_pt2);
  }

  void Line::check_and_raise(const Point_t& p1, const Point_t& p2) const
  { 
    if(p1.size()!=3) throw GeoAlgoException("<<check_and_raise>> Pt1 is not 3 dimensional point!");
    if(p2.size()!=3) throw GeoAlgoException("<<check_and_raise>> Pt2 is not 3 dimensional point!");
    if(p1 == p2) throw GeoAlgoException("<<check_and_raise>> Two identical points not allowed for Line ctor!"); 
  }

}
#endif
/** @} */ // end of doxygen group 

