#ifndef BASICTOOL_GEOHALFLINE_CXX
#define BASICTOOL_GEOHALFLINE_CXX

#include "GeoHalfLine.h"
namespace geoalgo {

  HalfLine::HalfLine()
    : _start(3)
    , _dir(3)
  {Normalize();}

  HalfLine::HalfLine(const double x,    const double y,    const double z,
		     const double dirx, const double diry, const double dirz)
    : _start (x,    y,    z   )
    , _dir   (dirx, diry, dirz)
  {Normalize();}

  HalfLine::HalfLine(const Point_t& start, const Vector_t& dir)
    : _start ( start )
    , _dir   ( dir   )
  { 
    if(start.size()!=3 || dir.size()!=3)
      throw GeoAlgoException("HalfLine ctor accepts only 3D Point!");
    Normalize();
  }

  const Point_t& HalfLine::Start() const { return _start; }
  
  const Vector_t& HalfLine::Dir() const { return _dir; }
  
  void HalfLine::Start(const double x, const double y, const double z)
  { _start[0] = x; _start[1] = y; _start[2] = z; }

  void HalfLine::Dir(const double x, const double y, const double z)
  { _dir[0] = x; _dir[1] = y; _dir[2] = z; Normalize(); }

  void HalfLine::Start(const TVector3& pt)
  { _start[0] = pt[0]; _start[1] = pt[1]; _start[2] = pt[2]; }

  void HalfLine::Dir(const TVector3& dir)
  { _dir[0] = dir[0]; _dir[1] = dir[1]; _dir[2] = dir[2]; Normalize(); }

  void HalfLine::Normalize()
  {
    auto l = _dir.Length();
    if(!l)
      throw GeoAlgoException("<<Normalize>> cannot normalize 0-length direction vector!");
    
    // inf check commented out till compatible solution found... --kazu
    //if(isnan(l))
    //throw GeoAlgoException("<<Normalize>> cannot normalize inf-length direction vector!");
    _dir /= l;
  }
}
#endif


