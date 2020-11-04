#ifndef BASICTOOL_GEOCYLINDER_CXX
#define BASICTOOL_GEOCYLINDER_CXX

#include "GeoCylinder.h"

namespace geoalgo {

  Cylinder::Cylinder() 
    : Line()
    , _radius (0.)
  {}
  
  Cylinder::Cylinder(const double x_min, const double y_min, const double z_min,
		     const double x_max, const double y_max, const double z_max,
		     const double radius)
    : Line(x_min, y_min, z_min, x_max, y_max, z_max)
    , _radius ( radius )
  {}
    
  Cylinder::Cylinder(const Point_t& min, const Vector_t& max, const double radius)
    : Line(min, max)
    , _radius ( radius )
  { 
    if(min.size()!=3 || max.size()!=3)
      throw GeoAlgoException("Cylinder ctor accepts only 3D Point!");
  }
  
  bool Cylinder::Contain(const Point_t &pt) const {

    // get a vector that defines the axis of the cylinder
    Vector_t axis = _pt1-_pt2;
    Vector_t dirpt = pt-_pt2;

    // angle of point w.r.t. the axis
    double angleMin = axis.Angle(dirpt);
    
    // if the angle is > 90 -> outside -> return
    if (angleMin > 0.5*3.14)
      return false;
    
    // revert the axis direction
    axis = _pt2-_pt1;
    dirpt = pt-_pt1;
    angleMin = axis.Angle(dirpt);
    
    // if the angle is > 90 -> outside -> return
    if (angleMin > 0.5*3.14)
      return false;

    // if still here, all that is left to verify is
    // that the point isn't more than a radius
    // away from the cylinder axis
    // 1) make a line corresponding to the axis
    // 2) get the distance between the point and the line
    double radial_dist_sq = _geoAlgo.SqDist(*this,pt);
    
    if (radial_dist_sq > _radius*_radius)
      return false;
    
    return true;

  }
}
#endif


