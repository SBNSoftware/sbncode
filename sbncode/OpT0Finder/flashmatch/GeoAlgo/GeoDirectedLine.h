/**
 * \file GeoDirectedLine.h
 *
 * \ingroup GeoAlgo
 * 
 * \brief Class def header for a class DirectedLine
 *
 * @author David Caratelli
 */

/** \addtogroup GeoAlgo

    @{*/
#ifndef BASICTOOL_GEODIRECTEDLINE_H
#define BASICTOOL_GEODIRECTEDLINE_H

#include "GeoLine.h"

namespace geoalgo {
  /**
     \class DirectedLine
     @brief Representation of a 3D infinite line.
     Defines an infinite 3D line with a point and a direction.
     Line points are constructed like this:
     (pt, dir) -> (pt, pt+dir)
     It hides the point attributes from users for protecting the dimensionality.
  */
  class DirectedLine : public Line {

  public:

    /// Default ctor
    DirectedLine();

    /// Alternative ctor (1)
    DirectedLine(const double x, const double y, const double z,
		 const double dirx, const double diry, const double dirz);
    
    /// Altenartive ctor (2)
    DirectedLine(const Point_t& pt, const Vector_t& dir);

    /// Alternative ctor (3)
    DirectedLine(const HalfLine& l);

    /// Alternative ctor using template (3)
    template <class T, class U> DirectedLine(const T& pt, const U& dir)
      : Line(Point_t(pt), Point_t(pt+dir))
    {}
    
    Vector_t Dir() const;
    
  };

  typedef DirectedLine DirectedLine_t;
}
#endif
/** @} */ // end of doxygen group 

