/**
 * \file GeoSphere.h
 *
 * \ingroup GeoAlgo
 * 
 * \brief Class def header for a class HalfLine
 *
 * @author kazuhiro
 */

/** \addtogroup GeoAlgo
    
    @{*/
#ifndef BASICTOOL_GEOSPHERE_H
#define BASICTOOL_GEOSPHERE_H

#include "GeoVector.h"

namespace geoalgo {
  /**
     \class Spehere
     @brief Representation of a 3D sphere
     Defines a 3D Sphere having an center (Point_t) and a radius (double) 
  */
  class Sphere {
    
  public:

    Sphere();           ///< Default ctor
    virtual ~Sphere(){} ///< Default dtor

    /// Alternative ctor (0)
    Sphere(const double& x, const double& y, const double& z, const double& r);

    /// Altenartive ctor (1) - 1 Point
    Sphere(const Point_t& center, const double r=0);
    
    /// Alternative ctor (2) - 2 Points
    Sphere(const Point_t& pt1, const Point_t& pt2);

    /// Alternative ctor (3) - 3 Points
    Sphere(const Point_t& A, const Point_t& B, const Point_t& C);
    
    //  Alternative ctor (4) - 4 Points
    Sphere(const Point_t& A, const Point_t& B, const Point_t& C, const Point_t& D);

    // Alternative ctor (5) - Set of points
    Sphere(const std::vector< ::geoalgo::Point_t>& pts);

    //
    // Getters
    //
    const Point_t& Center() const; ///< Center getter
    double Radius() const;   ///< Radius getter

    //
    // Setters
    //
    void Center(const double x, const double y, const double z); ///< Center setter
    void Center(const Point_t& pt); ///< Center setter
    void Radius(const double& r);   ///< Radius setter

    // 
    // Utilities
    //
    bool Contain(const Point_t& p) const; ///< Judge if a point is contained within a sphere

  protected:

    void compat(const Point_t& p, const double r=0) const; ///< 3D point compatibility check
    void compat(const double& r) const; ///< Positive radius compatibility check

    /// Center of Sphere
    Point_t  _center;

    /// Radius of Sphere
    double _radius;

  public:

    // 
    // Templates
    // 
    /*
#ifndef __CINT__ // Not sure why but CINT has a problem with this ctor. FIXME
    template <class T> Sphere(const T& center, const double r=0)
      : Sphere(Point_t(center),r)
    {}
#endif
    */
    template <class T> Sphere(const T& pt1, const T& pt2)
      : Sphere(Point_t(pt1), Point_t(pt2))
    {}

    template <class T> Sphere(const T& A, const T& B, const T& C)
      : Sphere(Point_t(A), Point_t(B), Point_t(C))
    {}

    template <class T> Sphere(const T& A, const T& B, const T& C, const T& D)
      : Sphere(Point_t(A), Point_t(B), Point_t(C), Point_t(D))
    {}
    
    template <class T> Sphere(const std::vector<T>& pts)
    {
      std::vector< ::geoalgo::Vector> geo_pts;
      geo_pts.reserve(pts);
      for(auto const& p : pts) geo_pts.emplace_back(p);
      (*this) = Sphere(geo_pts);
    }

    template <class T> void Center(const T& pt) 
    { Center(Point_t(pt)); }

    template <class T> bool Contain(const T& p) const
    { return Contain(Point_t(p)); }

  };
  
  typedef Sphere Sphere_t;
}
#endif
/** @} */ // end of doxygen group 

