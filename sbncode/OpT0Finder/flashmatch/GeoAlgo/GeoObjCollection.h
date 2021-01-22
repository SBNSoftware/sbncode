/**
 * \file GeoObjCollection.h
 *
 * \ingroup GeoAlgo
 * 
 * \brief Class def header for a class GeoObjCollection
 *
 * @author kazuhiro
 */

/** \addtogroup GeoAlgo

    @{*/
#ifndef BASICTOOL_GEOOBJCOLLECTION_H
#define BASICTOOL_GEOOBJCOLLECTION_H

#include <iostream>
#include "GeoTrajectory.h"
#include "GeoAABox.h"
#include "GeoLineSegment.h"
#include "GeoHalfLine.h"
#include "GeoCone.h"
#include "GeoSphere.h"
#include "GeoAlgoException.h"
#include <map>
namespace geoalgo {

  /**
     \class GeoObjCollection
  */
  class GeoObjCollection {
    
  public:
    
    /// Default constructor
    GeoObjCollection();
    
    /// Default destructor
    virtual ~GeoObjCollection(){}

    void Clear();
    
    void Add(const Point_t& pt, std::string label="", std::string c="");
	     
    void Add(const AABox_t& box, std::string label="", std::string c="");
	     
    void Add(const LineSegment_t& seg, std::string label="", std::string c="");

    void Add(const HalfLine_t& seg, std::string label="", std::string c="");
	     
    void Add(const Trajectory_t& trj, std::string label="", std::string c="");

    void Add(const Cone_t& cone, std::string label="", std::string c="");

    void Add(const Sphere_t& sphere, std::string label="", std::string c="");

    const std::vector< geoalgo::Point_t >& Point() const { return _pt_v; }
    const std::vector< std::string >& PointColor() const { return _pt_col; }

    const std::vector< geoalgo::AABox_t >& AABox() const { return _box_v; }
    const std::vector< std::string >& AABoxColor() const { return _box_col; }

    const std::vector< geoalgo::LineSegment_t >& LineSegment() const { return _seg_v; } 
    const std::vector< std::string >& LineSegmentColor() const { return _seg_col; } 

    const std::vector< geoalgo::HalfLine_t >& HalfLine() const { return _lin_v; } 
    const std::vector< std::string >& HalfLineColor() const { return _lin_col; } 

    const std::vector< geoalgo::Trajectory_t >& Trajectory() const { return _trj_v; }
    const std::vector< std::string >& TrajectoryColor() const { return _trj_col; }

    const std::vector< geoalgo::Cone_t >& Cone() const { return _cone_v; }
    const std::vector< std::string >& ConeColor() const { return _cone_col; }

    const std::vector< geoalgo::Sphere_t >& Sphere() const { return _sphere_v; }
    const std::vector< std::string >& SphereColor() const { return _sphere_col; }

    const std::map<geoalgo::Point_t,std::string>& Labels() const  { return _labels;}

  protected:

    const Point_t& _Point_(size_t i) const
    { return _pt_v[i]; }

    const AABox_t& _AABox_(size_t i) const
    { return _box_v[i]; }

    const LineSegment_t& _LineSegment_(size_t i) const
    { return _seg_v[i]; }

    const Trajectory_t& _Trajectory_(size_t i) const
    { return _trj_v[i]; }

    const Cone_t& _Cone_(size_t i) const
    { return _cone_v[i]; }

    const Sphere_t& _Sphere_(size_t i) const
    { return _sphere_v[i]; }

    void _AddLabel_(const Point_t& pt,
		    std::string label);

    std::vector< geoalgo::Point_t       > _pt_v;
    std::vector< std::string            > _pt_col; 
    std::vector< geoalgo::AABox_t       > _box_v;
    std::vector< std::string            > _box_col; 
    std::vector< geoalgo::LineSegment_t > _seg_v;
    std::vector< std::string            > _seg_col; 
    std::vector< geoalgo::HalfLine_t    > _lin_v;
    std::vector< std::string            > _lin_col; 
    std::vector< geoalgo::Trajectory_t  > _trj_v;
    std::vector< std::string            > _trj_col; 
    std::vector< geoalgo::Cone_t        > _cone_v;
    std::vector< std::string            > _cone_col; 
    std::vector< geoalgo::Sphere        > _sphere_v;
    std::vector< std::string            > _sphere_col; 
    std::map<geoalgo::Point_t,std::string > _labels;

  };
  
}
#endif
/** @} */ // end of doxygen group 

