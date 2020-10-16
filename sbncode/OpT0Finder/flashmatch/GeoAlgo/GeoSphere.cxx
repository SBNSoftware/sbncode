#ifndef BASICTOOL_GEOSPHERE_CXX
#define BASICTOOL_GEOSPHERE_CXX

#include "GeoSphere.h"

namespace geoalgo {

  Sphere::Sphere() : _center (3)
		   , _radius (0)
  { for(auto& v : _center) v=0; }
    
  Sphere::Sphere(const double& x,const double& y,const double& z,const double& r)
    : Sphere()
  { _center[0] = x; _center[1] = y; _center[2] = z; _radius = r; }

  Sphere::Sphere(const Point_t& center, const double r)
    : _radius (r)
  {_center = center;}

  Sphere::Sphere(const Point_t& pt1, const Point_t& pt2)
  {
    compat(pt1);
    compat(pt2);
    _center = (pt1+pt2)/2.;
    _radius = pt1.Dist(pt2)/2.;
  }

  //  3-Point Constructor
  //  Real-Time Collision Blog
  //  http://realtimecollisiondetection.net/blog/?p=20
  Sphere::Sphere(const Point_t& A, const Point_t& B, const Point_t& C)
  {
    compat(A);
    compat(B);
    compat(C);
    // any three points are co-planar
    // (if collinear no sphere passing  through all 3)
    // These 3 points make a triangle
    // find the perpendicular bi-sectors to the segments
    // making up the triangle. They will intersect
    // at the sphere's center
    
    // check if collinear. If so return exception
    Vector_t AB(B-A);
    Vector_t AC(C-A);
    Vector_t BC(C-B);
    
    double dABAB = AB.Dot(AB);
    double dACAC = AC.Dot(AC);
    double dABAC = AB.Dot(AC);
    
    double d = dABAB * dACAC - dABAC * dABAC;
    double s,t;
    
    // if d == 0 they lie on one line
    if (d == 0){
      std::cout << "d is 0!" << std::endl;
      double lenAB = AB.Length();
      double lenAC = AC.Length();
      double lenBC = BC.Length();
      // which segment is longest?
      if ( (lenAB > lenAC) && (lenAB > lenBC) ){
	_center = (A+B)/2.;
	_radius = _center.Dist(A);
      }
      else if( lenAC > lenBC ){
	_center = (A+C)/2.;
	_radius = _center.Dist(A);
      }
      else{
	_center = (B+C)/2;
	_radius = _center.Dist(B);
      }
    }// if d == 0
    
    else{
      s = 0.5 * ( dABAB * dACAC - dACAC * dABAC ) / d;
      t = 0.5 * ( dACAC * dABAB - dABAB * dABAC ) / d;
      
      // if s & t both > 0 && 1-s-t also > 0 then P = A + s*(B-A) + t*(C-A) is the center
      if ( (s > 0) && (t > 0) && ((1-s-t) > 0) ){
	_center = A+(B-A)*s+(C-A)*t;
	_radius = _center.Dist(A);
      }
      
      // otherwise only one will be negative. The side it belongs on will be
      // the longest side and will determine the side to take as diameter
      else if (s <= 0){
	// side AB is the one
	_center = (A+C)/2.;
	_radius = _center.Dist(A);
      }
      else if (t <= 0){
	// AC is the side
	_center = (A+B)/2.;
	_radius = _center.Dist(A);
      }
      else{
	_center = (B+C)/2;
	_radius = _center.Dist(B);
      }
    }// else (if d not equal to 0)
    
  }
  
  //  Alternative ctor (4) - 4 Points
  //  Real-Time Collision Blog
  //  http://realtimecollisiondetection.net/blog/?p=20
  /*
  Sphere::Sphere(const Point_t& A, const Point_t& B, const Point_t& C, const Point_t& D){

    compat(A);
    compat(B);
    compat(C);
    compat(D);

    // get sphere from 3 points (A,B,C)
    Vector_t AB(B-A);
    Vector_t AC(C-A);
    Vector_t AD(D-A);
    
    double dABAB = AB.Dot(AB);
    double dACAC = AC.Dot(AC);
    double dADAD = AD.Dot(AD);
    double dABAC = AB.Dot(AC);
    double dABAD = AB.Dot(AD);
    double dACAD = AC.Dot(AD);
    
    double d = 4*dABAC*dABAD*dACAD;

    if (d==0)
      throw GeoAlgoException("GeoSphere Exception: I think it means 3 points collinear. Find out which and call 3 point constructor - TO DO");
    
    double s = (dABAC*dACAD*dADAD + dABAD*dACAC*dACAD - dABAB*dACAD*dACAD)/d;
    double t = (dABAB*dACAD*dABAD + dABAD*dABAC*dADAD - dABAD*dABAD*dACAC)/d;
    double u = (dABAB*dABAC*dACAD + dABAC*dABAD*dACAC - dABAC*dABAC*dADAD)/d;
    
    // if everything positive! P = A + s(B-A) + t(C-A) + u(D-A)
    if ( (s > 0) && (t > 0) && (u > 0) && ((1-s-t-u) > 0) ){
      _center = A + AB*s + AC*t + AD*u;
      _radius = _center.Dist(A);
    }

    std::cout << "the center size at this stage is: " << _center.size() << std::endl;
    // TEMPORARY
    // otherwise find the 4 possible sphere combinations,
    // which contains the 4th point,
    // and if multiple ones choose the one with the smallest radius
    Sphere tmp = Sphere(A,B,C);
    _radius = kINVALID_DOUBLE;
    if (tmp.Contain(D)){
      _center = tmp.Center();
      _radius = tmp.Radius();
    }
    tmp = Sphere(A,B,D);
    if (tmp.Contain(C)){
      if (tmp.Radius() < _radius){
	_center = tmp.Center();
	_radius = tmp.Radius();
      }
    }
    tmp = Sphere(A,C,D);
    if (tmp.Contain(B)){
	if (tmp.Radius() < _radius){
	  _center = tmp.Center();
	  _radius = tmp.Radius();
	}
    }
    tmp = Sphere(B,C,D);
    if (tmp.Contain(A)){
      if (tmp.Radius() < _radius){
	_center = tmp.Center();
	_radius = tmp.Radius();
      }
    }

    std::cout << "the center size is: " << _center.size() << std::endl;

  }
  */

  // 4-point constructor
  Sphere::Sphere(const Point_t& A, const Point_t& B, const Point_t& C, const Point_t& D){

    compat(A);
    compat(B);
    compat(C);
    compat(D);

    // let's make sure there aren't duplicates...if so -> call a
    // different constructor
    std::vector<geoalgo::Point_t> valid_points = {A};
    bool duplicate = false;
    for (auto const& pt : valid_points){
      if (pt.SqDist(B) < 0.0001)
	duplicate = true;
    }
    if (duplicate == false)
      valid_points.push_back(B);
    duplicate = false;
    for (auto const& pt : valid_points){
      if (pt.SqDist(C) < 0.0001)
	duplicate = true;
    }
    if (duplicate == false)
      valid_points.push_back(C);
    duplicate = false;
    for (auto const& pt : valid_points){
      if (pt.SqDist(D) < 0.0001)
	duplicate = true;
    }
    if (duplicate == false)
      valid_points.push_back(D);
    
    // if we have less then 4 points -> call the appropriate constructor
    if (valid_points.size() < 4){
      (*this) = Sphere(valid_points);
      return;
    }

    // get sphere from 3 points (A,B,C)
    Vector_t AB(B-A);
    Vector_t AC(C-A);
    Vector_t AD(D-A);

    double dABAB = AB.Dot(AB);
    double dACAC = AC.Dot(AC);
    double dADAD = AD.Dot(AD);
    double dABAC = AB.Dot(AC);
    double dABAD = AB.Dot(AD);
    double dACAD = AC.Dot(AD);
    
    double d = 4*dABAC*dABAD*dACAD;

    if (d==0){
      // are any points duplicates? if so 
      // find the points that are collinear and call constructor
      // for the
      throw GeoAlgoException("GeoSphere Exception: I think it means 3 points collinear. Find out which and call 3 point constructor - TO DO");
    }
    
    double s = (dABAC*dACAD*dADAD + dABAD*dACAC*dACAD - dABAB*dACAD*dACAD)/d;
    double t = (dABAB*dACAD*dABAD + dABAD*dABAC*dADAD - dABAD*dABAD*dACAC)/d;
    double u = (dABAB*dABAC*dACAD + dABAC*dABAD*dACAC - dABAC*dABAC*dADAD)/d;
    
    // if everything positive! P = A + s(B-A) + t(C-A) + u(D-A)
    if ( (s > 0) && (t > 0) && (u > 0) && ((1-s-t-u) > 0) ){
      _center = A + AB*s + AC*t + AD*u;
      _radius = _center.Dist(A);
    }
    else{
      // take the largest side and use it as the diameter
      double maxdist = A.Dist(B);
      Vector_t max1 = A;
      Vector_t max2 = B;
      if (A.Dist(C) > maxdist){
	maxdist = A.Dist(C);
	max1 = A;
	max2 = C;
      }
      if (A.Dist(D) > maxdist){
	maxdist = A.Dist(D);
	max1 = A;
	max2 = D;
      }
      if (B.Dist(C) > maxdist){
	maxdist = B.Dist(C);
	max1 = B;
	max2 = C;
      }
      if (B.Dist(D) > maxdist){
	maxdist = B.Dist(D);
	max1 = B;
	max2 = D;
      }
      if (C.Dist(D) > maxdist){
	maxdist = C.Dist(D);
	max1 = C;
	max2 = D;
      }
      _center = (max1+max2)/2.;
      _radius = max1.Dist(max2)/2.;
    }
      

    // TEMPORARY
    // otherwise find the 4 possible sphere combinations,
    // which contains the 4th point,
    // and if multiple ones choose the one with the smallest radius
    Sphere tmp = Sphere(A,B,C);
    //_radius = kINVALID_DOUBLE;
    if (tmp.Contain(D)){
      _center = tmp.Center();
      _radius = tmp.Radius();
    }
    tmp = Sphere(A,B,D);
    if (tmp.Contain(C)){
      if (tmp.Radius() < _radius){
	_center = tmp.Center();
	_radius = tmp.Radius();
      }
    }
    tmp = Sphere(A,C,D);
    if (tmp.Contain(B)){
	if (tmp.Radius() < _radius){
	  _center = tmp.Center();
	  _radius = tmp.Radius();
	}
    }
    tmp = Sphere(B,C,D);
    if (tmp.Contain(A)){
      if (tmp.Radius() < _radius){
	_center = tmp.Center();
	_radius = tmp.Radius();
      }
    }

  }
  
  
  // Alternative ctor (5) - Set of points
  Sphere::Sphere(const std::vector< ::geoalgo::Point_t>& pts)
    : _center(0,0,0)
    , _radius(0)
  {

    switch(pts.size()) {
    case 0:
      break;
    case 1: _center = pts.front();
      break;
    case 2: (*this) = Sphere(pts[0],pts[1]);
      break;
    case 3: (*this) = Sphere(pts[0],pts[1],pts[2]);
      break;
    case 4: (*this) = Sphere(pts[0],pts[1],pts[2],pts[3]);
      break;
    default:
      throw GeoAlgoException("Cannot call Sphere constructor with more than 4 points. Something went wront");	
    }

  }
  
  const Point_t& Sphere::Center() const { return _center; }
  
  double Sphere::Radius() const { return _radius; }
  
  void Sphere::Center(const double x, const double y, const double z)
  { _center[0] = x; _center[1] = y; _center[2] = z; }

  void Sphere::Center(const Point_t& center)
  { compat(center); _center = center; }

  void Sphere::Radius(const double& r)
  { compat(r); _radius = r; }

  bool Sphere::Contain(const Point_t& p) const
  {
    _center.compat(p);
    return ( p._Dist_(_center) < _radius );
  }
  
  void Sphere::compat(const Point_t& p, const double r) const
  { 
    if(p.size()!=3) throw GeoAlgoException("Only 3D points allowed for sphere"); 
    compat(r);
  }
  
  void Sphere::compat(const double& r) const
  { if(r<0) throw GeoAlgoException("Only positive value allowed for radius"); }
  
}

#endif


