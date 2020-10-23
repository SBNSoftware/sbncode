/**
 * \file GeoCone.h
 *
 * \ingroup GeoAlgo
 * 
 * \brief Class def header for a class HalfLine
 *
 * @author David Caratelli
 */

/** \addtogroup GeoAlgo

    @{*/
#ifndef BASICTOOL_GEOCONE_H
#define BASICTOOL_GEOCONE_H

#include "GeoVector.h"
#include "GeoHalfLine.h"
namespace geoalgo {
  /**
     \class Cone
     @brief Representation of a 3D semi-infinite line.
     Defines a 3D cone with the following properties:                                           \n
     Start point (or vertex), Direction, Length (or Length), Radius, opening angle              \n
     When 2 of Length, Radius, opening angle are defined the third is automatically set         
  */
  class Cone : public HalfLine {
    
  public:
    
    /// Default constructor
    Cone();

    /// Default destructor
    virtual ~Cone(){};

    /// Alternative ctor (1)
    Cone(const double x,    const double y,    const double z,
	 const double dirx, const double diry, const double dirz,
	 const double length, const double radius);

    /// Alternative ctor (2)
    Cone(const Point_t& start, const Vector_t& dir,
	 const double length, const double radius);

    //
    // Getters
    //
    double Length() const; ///< Length getter
    double Radius() const; ///< Length getter
    double Angle () const; ///< Angle getter

    //
    // Setters
    //
    void Length(const double l); ///< Length setter
    void Radius(const double r); ///< Radius setter

  protected:

    double _length; ///< Helight (length) of the cone
    double _radius; ///< Radius of the cone at the base
    double _angle;  ///< Opening Angle

  public:
    //
    // Template
    //
    /// Alternative ctor using template (3)
    template <class T, class U> Cone(const T& start, const U& dir)
      : Cone(Point_t(start), Vector_t(dir))
    {}

  };
  
  typedef Cone Cone_t;
}
#endif
/** @} */ // end of doxygen group 

