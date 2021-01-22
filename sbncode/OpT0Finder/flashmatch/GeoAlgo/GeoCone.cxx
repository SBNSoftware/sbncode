#ifndef BASICTOOL_GEOCONE_CXX
#define BASICTOOL_GEOCONE_CXX

#include "GeoCone.h"
#include <sstream>
namespace geoalgo {

  Cone::Cone() : HalfLine()
  {
    _length = 1;
    _radius = 1;
    _angle = atan(_radius/_length);
  }
    
  Cone::Cone(const double x,    const double y,    const double z,
	     const double dirx, const double diry, const double dirz,
	     const double length, const double radius)
    : HalfLine(x, y, z, dirx, diry, dirz)
  {
    if (length == 0){
      std::ostringstream msg;
      msg << "<<" << __FUNCTION__ << ">>"
	  << " Cone Length cannot be 0." << std::endl;
      throw GeoAlgoException(msg.str());
    }
    _length = length;
    _radius = radius;
    _angle = atan(_radius/_length);
  }
  
  Cone::Cone(const Point_t& start, const Vector_t& dir,
	     const double length, const double radius)
    : HalfLine( start,  dir )
  { 
    if (length == 0){
      std::ostringstream msg;
      msg << "<<" << __FUNCTION__ << ">>"
	  << " Cone Length cannot be 0." << std::endl;
      throw GeoAlgoException(msg.str());
    }
    _length = length;
    _radius = radius;
    _angle  = atan(_radius/_length);
  }

  double Cone::Length() const { return _length; }

  double Cone::Radius() const { return _radius; }

  double Cone::Angle() const { return _angle; }
    
  void Cone::Length(const double l) 
  { 
    if (l == 0){
      std::ostringstream msg;
      msg << "<<" << __FUNCTION__ << ">>"
	  << " Cone Length cannot be 0." << std::endl;
      throw GeoAlgoException(msg.str());
    }
    _length = l; 
    _angle = atan(_radius/_length);
  }
  
  void Cone::Radius(const double r) { _radius = r; _angle = atan(_radius/_length); }
}
#endif


