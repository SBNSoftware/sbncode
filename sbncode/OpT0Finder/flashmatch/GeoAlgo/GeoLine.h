/**
 * \file GeoLine.h
 *
 * \ingroup GeoAlgo
 * 
 * \brief Class def header for a class Line
 *
 * @author David Caratelli
 */

/** \addtogroup GeoAlgo

    @{*/
#ifndef BASICTOOL_GEOLINE_H
#define BASICTOOL_GEOLINE_H

#include "GeoHalfLine.h"

namespace geoalgo {
  /**
     \class Line
     @brief Representation of a 3D infinite line.
     Defines an infinite 3D line by having 2 points which completely determine the line
     along which the line extends. It hides the point attributes from users for   \n
     protecting the dimensionality.
  */
  class Line {
    
  public:
    
    /// Default constructor
    Line();

    /// Default destructor
    virtual ~Line(){}

    /// Alternative ctor (1)
    Line(const double x1, const double y1, const double z1,
	 const double x2, const double y2, const double z2);

    /// Altenartive ctor (2)
    Line(const Point_t& pt1, const Point_t& pt2);

    //
    // Getters
    //
    const Point_t& Pt1() const; ///< Start getter
    const Point_t& Pt2() const; ///< Direction getter

    //
    // Setters
    //
    void Pt1(const double x, const double y, const double z); ///< Pt1 setter
    void Pt2(const double x, const double y, const double z); ///< Pt2 setter

  protected:

    /// Compatibility check
    void check_and_raise(const Point_t& p1, const Point_t& p2) const;

    Point_t  _pt1; ///< First point denoting infinite line
    Vector_t _pt2; ///< Second point denoting infinite line

  public:
    //
    // Template
    //
    /// Alternative ctor using template (3)
    template <class T, class U> Line(const T& pt1, const U& pt2)
      : Line(Point_t(pt1), Point_t(pt2))
    {}

    /// Pt1 setter template
    template<class T>
    void Pt1(const T& pt1)
    { 
      _pt1 = Point_t(pt1); 
      check_and_raise(_pt1,_pt2);
    }
    
    /// Pt2 setter template
    template<class T>
    void Pt2(const T& pt2)
    { 
      _pt2 = Vector_t(pt2);
      check_and_raise(_pt1,_pt2);
    }
  };


  typedef Line Line_t;
}
#endif
/** @} */ // end of doxygen group 

