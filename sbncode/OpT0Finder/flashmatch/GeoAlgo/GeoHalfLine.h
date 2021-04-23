/**
 * \file GeoHalfLine.h
 *
 * \ingroup GeoAlgo
 * 
 * \brief Class def header for a class HalfLine
 *
 * @author kazuhiro
 */

/** \addtogroup GeoAlgo

    @{*/
#ifndef BASICTOOL_GEOHALFLINE_H
#define BASICTOOL_GEOHALFLINE_H

#include "GeoVector.h"
namespace geoalgo {
  /**
     \class HalfLine
     @brief Representation of a 3D semi-infinite line.
     Defines a semi-infinite 3D line by having a start point (Point_t) and direction (Vector_t) \n
     along which the line extends. It hides the start and direction attributes from users for   \n
     protecting the dimensionality.
  */
  class HalfLine {
    
  public:
    
    /// Default constructor
    HalfLine();

    /// Default destructor
    virtual ~HalfLine(){};

    /// Alternative ctor (1)
    HalfLine(const double x,    const double y,    const double z,
	     const double dirx, const double diry, const double dirz);

    /// Altenartive ctor (2)
    HalfLine(const Point_t& start, const Vector_t& dir);
    
    const Point_t&  Start () const; ///< Start getter
    const Vector_t& Dir   () const; ///< Direction getter
    
    void Start(const double x, const double y, const double z); ///< Start setter
    void Dir  (const double x, const double y, const double z); ///< Dir setter

    void Start(const TVector3& pt ); ///< Start setter
    void Dir  (const TVector3& dir); ///< Dir setter

  protected:

    void Normalize(); ///< Normalize direction
    Point_t  _start;  ///< Beginning of the half line
    Vector_t _dir;    ///< Direction of the half line from _start

  public:

    //
    // Template
    // 

    /// Alternative ctor using template (3)
    template <class T, class U> HalfLine(const T& start, const U& dir)
      : HalfLine(Point_t(start), Vector_t(dir))
    {}

    /// Start setter template
    template<class T>
    void Start(const T& pos)
    { 
      _start = Point_t(pos); 
      if(_start.size()!=3) throw GeoAlgoException("<<Start>> Only 3 dimensional start point allowed!"); 
    }
    
    /// Dir setter template
    template<class T>
    void Dir(const T& dir)
    { 
      _dir = Vector_t(dir);
      if(_dir.size()!=3) throw GeoAlgoException("<<Start>> Only 3 dimensional start point allowed!"); 
      Normalize();
    }

  };
  
  typedef HalfLine HalfLine_t;
}
#endif
/** @} */ // end of doxygen group 

