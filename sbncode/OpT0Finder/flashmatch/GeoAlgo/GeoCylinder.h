/**
 * \file GeoCylinder.h
 *
 * \ingroup GeoAlgo
 * 
 * \brief Class def header for a class Cylinder
 *
 * @author david caratelli
 */

/** \addtogroup GeoAlgo

    @{*/
#ifndef BASICTOOL_GEOCYLINDER_H
#define BASICTOOL_GEOCYLINDER_H

#include "GeoAlgo.h"
#include "GeoLine.h"

namespace geoalgo {
  /**
     \class Cylinder
     @brief Representation of a 3D Cylinder volume.
     A Cylinder object inherits from a geoalgo::Line
     @input 2 points, which define the line representing
     the central axis of the cylinder
     @input a radius, defining the radius of the cylinder
  */
  class Cylinder : public Line {
    
  public:
    
    /// Default constructor
    Cylinder();

    /// Default destructor
    virtual ~Cylinder(){};
    
    /// Alternative ctor (0)
    Cylinder(const double x_min, const double y_min, const double z_min,
	     const double x_max, const double y_max, const double z_max,
	     const double radius);
    
    /// Altenartive ctor (1)
    Cylinder(const Point_t& min, const Vector_t& max, const double radius);
    
    /// Containment evaluation
    bool Contain(const Point_t &pt) const; ///< Test if a point is contained within the box

    /// Getters
    double GetRadius() { return _radius; }
    /// Setters
    void SetRadius(double r) { _radius = r; }
    
  protected:
    
    double  _radius; ///< Radius of the cylinder

    // geoalgo utility
    GeoAlgo _geoAlgo;

  };

  typedef Cylinder Cylinder_t;
}
#endif
/** @} */ // end of doxygen group 

