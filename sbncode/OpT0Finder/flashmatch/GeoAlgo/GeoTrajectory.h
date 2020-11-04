/**
 * \file GeoObjects.h
 *
 * \ingroup GeoAlgo
 * 
 * \brief Class def header for a class Trajectory
 *
 * @author kazuhiro
 */

/** \addtogroup GeoAlgo

    @{*/
#ifndef BASICTOOL_GEOTRAJECTORY_H
#define BASICTOOL_GEOTRAJECTORY_H

#include "GeoVector.h"

namespace geoalgo {

  /**
     \class Trajectory
     This class represents a trajectory which is an ordered list of Point.
     It is a friend class w/ geoalgo::Point_t hence it has an access to protected functions that avoids
     dimensionality sanity checks for speed.
   */
  class Trajectory : public std::vector<geoalgo::Vector> {

  public:
    
    /// Default ctor to specify # points and dimension of each point
    Trajectory(size_t npoints=0, size_t ndimension=0);

    /// Default dtor
    virtual ~Trajectory(){}

    /// Alternative ctor (0) using a vector of mere vector point expression
    Trajectory(const std::vector<std::vector<double> > &obj);

    /// Alternative ctor (1) using a vector of point
    Trajectory(const std::vector<geoalgo::Point_t> &obj);

    //
    // Getters
    //
    double Length(size_t start_step=0,size_t end_step=0) const; ///< The summed-length along all trajectory points
    bool   IsLonger(double) const; ///< Check if the trajectory is longer than specified value
    Vector Dir(size_t i=0)  const; ///< The direction at a specified trajectory point

    //
    // Setters
    //
    void push_back(const Point_t& obj); ///< push_back overrie w/ dimensionality check 

    inline Trajectory& operator+=(const Point_t& rhs)
    { push_back(rhs); return *this; }

    //
    // utility
    //
    void compat(const Point_t& obj)    const; ///< Dimensionality check function w/ Trajectory
    void compat(const Trajectory &obj) const; ///< Dimensionality check function w/ Point_t

  protected:

    /// Returns a direction vector at a specified trajectory point w/o size check
    Vector _Dir_(size_t i) const;

  public:

    //
    // templates
    //
    /// push_back template
    template <class T>
    void push_back(const T& obj) 
    { Point_t pt(obj); push_back(pt); }

  public:

    /// Streamer
#ifndef __CINT__
    friend std::ostream& operator << (std::ostream &o, Trajectory const& a)
    { o << "Trajectory with " << a.size() << " points " << std::endl;
      for(auto const& p : a )
	o << " " << p << std::endl;
      return o;
    }
#endif

  };

  typedef Trajectory Trajectory_t;

}

#endif
/** @} */ // end of doxygen group 

