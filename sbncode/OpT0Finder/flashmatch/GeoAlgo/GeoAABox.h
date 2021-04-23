/**
 * \file GeoAABox.h
 *
 * \ingroup GeoAlgo
 * 
 * \brief Class def header for a class AABox
 *
 * @author kazuhiro
 */

/** \addtogroup GeoAlgo

    @{*/
#ifndef BASICTOOL_GEOAABOX_H
#define BASICTOOL_GEOAABOX_H

#include "GeoHalfLine.h"

namespace geoalgo {
  /**
     \class AABox
     @brief Representation of a 3D rectangular box which sides are aligned w/ coordinate axis.
     A representation of an Axis-Aligned-Boundary-Box, a simple & popular representation of   \n
     3D boundary box for collision detection. The concept was taken from the reference,       \n
     Real-Time-Collision-Detection (RTCD), and in particular Ch. 4.2 (page 77):               \n

     Ref: http://realtimecollisiondetection.net

     This class uses one of the simplest representation for AABox: "min-max" representation.   \n
     Though this method requires storing 6 floating point values, class attributes (i.e.      \n
     "min" and "max" points) store intuitive values for most UB analyzers. Also it simplifies \n
     utility function implementations.
  */
  class AABox {
    
  public:
    
    /// Default constructor
    AABox();

    /// Default destructor
    virtual ~AABox(){};
    
    /// Alternative ctor (0)
    AABox(const double x_min, const double y_min, const double z_min,
	  const double x_max, const double y_max, const double z_max);
    
    /// Altenartive ctor (1)
    AABox(const Point_t& min, const Vector_t& max);

    //
    // Attribute accessor
    //
    const Point_t& Min() const; ///< Minimum point getter
    const Point_t& Max() const; ///< Maximum point getter
    void Min(const double x, const double y, const double z); ///< Minimum point setter
    void Max(const double x, const double y, const double z); ///< Maximum point setter
    bool Contain(const Point_t &pt) const; ///< Test if a point is contained within the box
    
  protected:
    
    Point_t _min; ///< Minimum point
    Point_t _max; ///< Maximum point

  public:

    //
    // Template
    //
    /// Alternative ctor using template (3)
    template <class T, class U> AABox(const T& min, const U& max)
      : AABox(Point_t(min), Point_t(max))
    {}

    /// Streamer
    #ifndef __CINT__
    friend std::ostream& operator << (std::ostream &o, ::geoalgo::AABox const& a)
    { o << "AABox Min " << a.Min() << " Max " << a.Max(); return o; }
    #endif
    
  };

  typedef AABox AABox_t;

}
#endif
/** @} */ // end of doxygen group 

