/**
 * \file GeoAlgo.h
 *
 * \ingroup GeoAlgo
 * 
 * \brief Class def header for a class GeoAlgo
 *
 * @author kazuhiro
 */

/** \addtogroup GeoAlgo

    @{*/
#ifndef BASICTOOL_GEOALGO_H
#define BASICTOOL_GEOALGO_H

#include <algorithm> // used for std::find
#include "GeoLine.h"
#include "GeoHalfLine.h"
#include "GeoLineSegment.h"
#include "GeoTrajectory.h"
#include "GeoCone.h"
#include "GeoAABox.h"
#include "GeoSphere.h"

namespace geoalgo {

  /**
     \class GeoAlgo
     @brief Algorithm to compute various geometrical relation among geometrical objects.
     In particular functions to inspect following relations are implemented: \n
     0) Distance between geometrical objects \n
     1) Closest point of approach            \n
     2) Intersection points                  \n
     3) Containment/Overlap of objects       \n
     4) Common Origin functions              \n
     5) Bounding Sphere functions            \n

     Most functions are taken from the reference Real-Time-Collision-Detection (RTCD):
     Ref: http://realtimecollisiondetection.net
  */
  class GeoAlgo{
    
  public:
    
    /// Default constructor
    GeoAlgo(){}
    
    /// Default destructor
    virtual ~GeoAlgo(){}

    //
    // Intersections
    //

    /// Intersection between a HalfLine and an AABox
    std::vector<Point_t> Intersection(const AABox_t& box, const HalfLine_t& line, bool back=false) const;
    /// Intersection between a HalfLine and an AABox
    std::vector<Point_t> Intersection(const HalfLine_t& line, const AABox_t& box, bool back=false) const
    { return Intersection(box, line, back); }

    /// Intersection between LineSegment and an AABox
    std::vector<Point_t> Intersection(const AABox_t& box, const LineSegment_t& l) const;
    /// Intersection between LineSegment and an AABox
    std::vector<Point_t> Intersection(const LineSegment_t& l, const AABox_t& box) const
    { return Intersection(box,l); }

    /// Intersection between Trajectory and an AABox
    std::vector<Point_t> Intersection(const AABox_t& box, const Trajectory_t& trj) const;
    /// Intersection between Trajectory and an AABox
    std::vector<Point_t> Intersection(const Trajectory_t& trj, const AABox_t& box) const
    { return Intersection(box,trj); }

    /// LineSegment sub-segment of HalfLine inside an AABox
    LineSegment_t BoxOverlap(const AABox_t& box, const HalfLine_t& line) const;
    /// LineSegment sub-segment of HalfLine inside an AABox
    LineSegment_t BoxOverlap(const HalfLine_t& line, const AABox_t& box) const
    { return BoxOverlap(box, line); }


    /// Get Trajectory inside box given some input trajectory -> now assumes trajectory cannot exit and re-enter box
    Trajectory_t BoxOverlap(const AABox_t& box, const Trajectory_t& trj) const;
    /// Get Trajectory inside box given some input trajectory -> now assumes trajectory cannot exit and re-enter box
    Trajectory_t BoxOverlap(const Trajectory_t& trj, const AABox_t& box) const
    { return BoxOverlap(box, trj); }
        

    //************************************************
    //CLOSEST APPROACH BETWEEN POINT AND INFINITE LINE
    //************************************************
    // min distance between point and infinite line
    double SqDist(const Line_t& line, const Point_t& pt) const
    { pt.compat(line.Pt1()); return _SqDist_(line,pt); }
    // min distance between point and infinite line
    double SqDist(const Point_t& pt, const Line_t& line) const
    { return SqDist(line,pt); }
    // closest point (on infinite line) from point
    Point_t ClosestPt(const Line_t& line, const Point_t& pt) const
    { pt.compat(line.Pt1()); return _ClosestPt_(pt,line); }
    // closest point (on infinite line) from point
    Point_t ClosestPt(const Point_t& pt, const Line_t& line) const
    { return _ClosestPt_(pt,line); }


    //*******************************************
    //CLOSEST APPROACH BETWEEN TWO INFINITE LINES
    //*******************************************
    // Closest approach between two infinite line segments - keep track of closest approach points
    double SqDist(const Line_t& l1, const Line_t& l2, Point_t& L1, Point_t& L2) const
    { l1.Pt1().compat(l2.Pt1()); return _SqDist_(l1, l2, L1, L2); }
    // Closest approach between two infinite line segments - don't keep track of closest approach points
    double SqDist(const Line_t& l1, const Line_t& l2) const
    { Point_t L1; Point_t L2; return SqDist(l1, l2, L1, L2); }
    

    //************************************************
    //CLOSEST APPROACH BETWEEN TWO HALF-INFINITE LINES
    //************************************************
    // Closest approach between two infinite line segments - keep track of closest approach points
    double SqDist(const HalfLine_t& l1, const HalfLine_t& l2, Point_t& L1, Point_t& L2) const
    { l1.Start().compat(l2.Start()); return _SqDist_(l1, l2, L1, L2); }
    // Closest approach between two infinite line segments - don't keep track of closest approach points
    double SqDist(const HalfLine_t& l1, const HalfLine_t& l2) const
    { Point_t L1; Point_t L2; return SqDist(l1, l2, L1, L2); }


    //******************************************
    //CLOSEST APPROACH BETWEEN TWO LINE SEGMENTS
    //******************************************
    /// LineSegment_t & LineSegment_t distance - keep track of points
    double SqDist(const LineSegment_t& seg1, const LineSegment_t& seg2, Point_t& c1, Point_t& c2) const
    { seg1.Start().compat(seg2.Start()); return _SqDist_(seg1, seg2, c1, c2); }
    /// LineSegment & LineSegment, don't keep track of points
    double SqDist(const LineSegment_t& seg1, const LineSegment_t& seg2) const
    { Point_t c1; Point_t c2; return SqDist(seg1, seg2, c1, c2); }


    //******************************************
    //CLOSEST APPROACH BETWEEN SEGMENT AND TRACK
    //******************************************
    /// LineSegment & Trajectory, keep track of points
    double SqDist(const LineSegment_t& seg, const Trajectory_t& trj, Point_t& c1, Point_t& c2) const;
    /// LineSegment & Trajectory, keep track of points
    double SqDist(const Trajectory_t& trj, const LineSegment_t& seg, Point_t& c1, Point_t& c2) const
    { return SqDist(seg, trj, c1, c2); }
    /// LineSegment & Trajectory, don't keep track of points
    double SqDist(const Trajectory_t& trj, const LineSegment_t& seg) const
    { Point_t c1; Point_t c2; return SqDist(seg, trj, c1, c2); }
    /// LineSegment & Trajectory, don't keep track of points
    double SqDist(const LineSegment_t& seg, const Trajectory_t& trj) const
    { Point_t c1; Point_t c2; return SqDist(seg, trj, c1, c2); }


    //****************************************
    //CLOSEST APPROACH BETWEEN TRACK AND TRACK
    //****************************************
    /// Trajectory & Trajectory, keep track of points
    double SqDist(const Trajectory_t& trj1, const Trajectory_t& trj2, Point_t& c1, Point_t& c2) const;
    /// Trajectory & Trajectory, don't keep track of points
    double SqDist(const Trajectory_t& trj1, const Trajectory_t& trj2) const
    { Point_t c1; Point_t c2; return SqDist(trj1, trj2, c1, c2); }


    //*****************************************************
    //CLOSEST APPROACH BETWEEN SEGMENT AND VECTOR OF TRACKS
    //*****************************************************
    /// LineSegment & vector of Trajectories, keep track of points
    double SqDist(const LineSegment_t& seg, const std::vector<geoalgo::Trajectory_t> &trj, Point_t& c1, Point_t& c2, int& trackIdx) const;
    /// LineSegment & vector of Trajectories, keep track of points
    double SqDist(const std::vector<geoalgo::Trajectory_t> &trj, const LineSegment_t& seg, Point_t& c1, Point_t& c2, int& trackIdx) const
    { return SqDist(seg, trj, c2, c1, trackIdx); }
    /// LineSegment & vector of Trajectories, don't keep track of points
    double SqDist(const std::vector<geoalgo::Trajectory_t> &trj, const LineSegment_t& seg) const
    { Point_t c1; Point_t c2; int trackIdx; return SqDist(seg, trj, c1, c2, trackIdx); }
    /// LineSegment & vector of Trajectories, don't keep track of points
    double SqDist(const LineSegment_t& seg, const std::vector<geoalgo::Trajectory_t> &trj) const
    { Point_t c1; Point_t c2; int trackIdx; return SqDist(seg, trj, c1, c2, trackIdx); }


    //*******************************************
    //CLOSEST APPROACH BETWEEN HALFLINE AND TRACK
    //*******************************************
    /// HalfLine & Trajectory, keep track of points
    double SqDist(const HalfLine_t& hline, const Trajectory_t& trj, Point_t& c1, Point_t& c2) const;
    /// HalfLine & Trajectory, keep track of points
    double SqDist(const Trajectory_t& trj, const HalfLine_t& hline, Point_t& c1, Point_t& c2) const
    { return SqDist(hline, trj, c2, c1); }
    /// HalfLine & Trajectory, don't keep track of points
    double SqDist(const Trajectory_t& trj, const HalfLine_t& hline) const
    { Point_t c1; Point_t c2; return SqDist(hline, trj, c1, c2); }
    /// HalfLine & Trajectory, don't keep track of points
    double SqDist(const HalfLine_t& hline, const Trajectory_t& trj) const
    { Point_t c1; Point_t c2; return SqDist(hline, trj, c1, c2); }


    //****************************************
    //CLOSEST APPROACH BETWEEN POINT AND TRACK
    //****************************************
    /// Point_t & Trajectory_t distance
    double SqDist(const Point_t& pt, const Trajectory_t& trj) const;
    /// Point_t & Trajectory_t distance
    double SqDist(const Trajectory_t& trj, const Point_t& pt) const
    { return SqDist(pt,trj); }
    /// Point_t & Trajectory_t closest point
    Point_t ClosestPt(const Point_t& pt, const Trajectory_t& trj) const
    { int idx=0; return ClosestPt(pt,trj,idx); }
    /// Point_t & Trajectory_t closest point
    Point_t ClosestPt(const Trajectory_t& trj, const Point_t& pt) const
    { int idx=0; return ClosestPt(pt,trj,idx); }
    /// Point_t & Trajectory_t closest point. Keep track of index of segment
    Point_t ClosestPt(const Point_t& pt, const Trajectory_t& trj, int& idx) const;
    /// Point_t & Trajectory_t closest point. Keep track of index of segment
    Point_t ClosestPt(const Trajectory_t& trj, const Point_t& pt, int& idx) const
    { return ClosestPt(pt,trj,idx); }


    //***************************************************
    //CLOSEST APPROACH BETWEEN POINT AND VECTOR OF TRACKS
    //***************************************************
    /// Point_t & Trajectory_t distance - keep track of which track
    double SqDist(const Point_t& pt, const std::vector<geoalgo::Trajectory_t> &trj, int& trackIdx) const;
    /// Point_t & Trajectory_t distance - keep track of which track
    double SqDist(const std::vector<geoalgo::Trajectory_t> &trj, const Point_t& pt, int& trackIdx) const
    { return SqDist(pt,trj, trackIdx); }
    /// Point_t & Trajectory_t distance - don't keep track
    double SqDist(const Point_t& pt, const std::vector<geoalgo::Trajectory_t> &trj) const
    { int trackIdx; return SqDist(pt, trj, trackIdx); }
    /// Point_t & Trajectory_t distance - don't keep track
    double SqDist(const std::vector<geoalgo::Trajectory_t> &trj, const Point_t& pt) const
    { int trackIdx; return SqDist(pt, trj, trackIdx); }
    /// Point_t & Trajectory_t closest point - keep track of which track is closest 
    Point_t ClosestPt(const Point_t& pt, const std::vector<geoalgo::Trajectory_t> &trj, int &trackIdx) const;
    /// Point_t & Trajectory_t closest point - keep track of which track is closest
    Point_t ClosestPt(const std::vector<geoalgo::Trajectory_t> &trj, const Point_t& pt, int &trackIdx) const
    { return ClosestPt(pt, trj, trackIdx); }
    /// Point_t & Trajectory_t closest point - don't keep track of which track is closest 
    Point_t ClosestPt(const Point_t& pt, const std::vector<geoalgo::Trajectory_t> &trj) const
    { int trackIdx; return ClosestPt(pt, trj, trackIdx); }
    /// Point_t & Trajectory_t closest point - don't keep track of which track is closest
    Point_t ClosestPt(const std::vector<geoalgo::Trajectory_t> &trj, const Point_t& pt) const
    { int trackIdx; return ClosestPt(pt, trj, trackIdx); }


    //***********************************************
    //CLOSEST APPROACH BETWEEN POINT AND LINE SEGMENT
    //***********************************************
    /// Point & LineSegment_t distance
    double SqDist(const Point_t& pt, const LineSegment_t& line) const
    { pt.compat(line.Start()); return _SqDist_(pt,line); }
    /// Point & LineSegment distance
    double SqDist(const LineSegment_t& line, const Point_t& pt) const
    { return SqDist(pt,line); }
    /// Point & LineSegment closest point
    Point_t ClosestPt(const Point_t& pt, const LineSegment_t& line) const
    { pt.compat(line.Start()); return _ClosestPt_(pt,line); }
    /// Point & LineSegment closest point
    Point_t ClosestPt(const LineSegment_t& line, const Point_t& pt) const
    { return ClosestPt(pt,line); }
    
    //********************************************
    //CLOSEST APPROACH BETWEEN POINT AND HALF LINE
    //********************************************
    /// Point & HalfLine distance
    double SqDist(const Point_t& pt, const HalfLine_t& line) const
    { pt.compat(line.Start()); return _SqDist_(pt,line); }
    /// Point & HalfLine distance
    double SqDist(const HalfLine_t& line, const Point_t& pt) const
    { return SqDist(pt,line); }
    /// Point & HalfLine closest point
    Point_t ClosestPt(const Point_t& pt, const HalfLine_t& line) const
    { pt.compat(line.Start()); return _ClosestPt_(pt,line); }
    /// Point & HalfLine closest point
    Point_t ClosestPt(const HalfLine_t& line, const Point_t& pt) const
    { return ClosestPt(pt,line); }


    //***************************************************
    //CLOSEST APPROACH BETWEEN HALF LINE AND LINE SEGMENT
    //***************************************************
    // half-line and line-segment. keep track of closest approach points
    double SqDist(const HalfLine_t& hline, const LineSegment_t& seg, Point_t& L1, Point_t& L2) const
    { hline.Start().compat(seg.Start()); return _SqDist_(hline, seg, L1, L2); }
    // half-line and line-segment. keep track of closest approach points
    double SqDist(const LineSegment_t& seg, const HalfLine_t& hline, Point_t& L1, Point_t& L2) const
    { return SqDist(hline,seg, L2, L1); }  
    // half-line and line-segment. Do not keep track of closest approach points
    double SqDist(const HalfLine_t& hline, const LineSegment_t& seg) const
    { Point_t L1; Point_t L2; return SqDist(hline, seg, L1, L2); }
    // half-line and line-segment. Do not keep track of closest approach points
    double SqDist(const LineSegment_t& seg, const HalfLine_t& hline) const
    { return SqDist(hline,seg); }  

    //***************************************************
    //CLOSEST APPROACH BETWEEN POINT AND AXIS ALIGNED BOX
    //***************************************************
    /// Point & AABox distance
    double SqDist(const Point_t& pt, const AABox_t& box) const
    { pt.compat(box.Min()); return _SqDist_(pt,box); }
    /// Point & AABox distance
    double SqDist(const AABox_t& box, const Point_t& pt)
    { return SqDist(pt,box); }
    /// Point & AABox closest point
    Point_t ClosestPt(const Point_t& pt, const AABox_t& box) const
    { pt.compat(box.Min()); return _ClosestPt_(pt,box); }
    /// Point & AABox closest point
    Point_t ClosestPt(const AABox_t& box, const Point_t& pt) const
    { return ClosestPt(pt,box); }


    //***************************************************************************
    //COMMON ORIGIN ALGORITHMS: DETERMINE IF TWO GEO-OBJECTS HAVE A COMMON ORIGIN
    //***************************************************************************
    /// Common origin: Line Segment & Line Segment. Do not keep track of origin
    double commonOrigin(const Line_t& lin1, const Line_t& lin2) const
    { Point_t origin(lin1.Pt1().size()); return commonOrigin(lin1,lin2, origin); }
    /// Common origin: Line Segment & Line Segment. Keep track of origin
    double commonOrigin(const Line_t& lin1, const Line_t& lin2, Point_t& origin) const
    { lin1.Pt1().compat(lin2.Pt1()); return _commonOrigin_(lin1, lin2, origin); }
    /// Common origin: Line Segment & Line Segment. Do not keep track of origin
    double commonOrigin(const LineSegment_t& seg1, const LineSegment_t& seg2, bool backwards=false) const
    { Point_t origin(seg1.Start().size()); return commonOrigin(seg1,seg2, origin, backwards); }
    /// Common origin: Line Segment & Line Segment. Keep track of origin
    double commonOrigin(const LineSegment_t& seg1, const LineSegment_t& seg2, Point_t& origin, bool backwards=false) const
    { seg1.Start().compat(seg2.Start()); return _commonOrigin_(seg1, seg2, origin, backwards); }
    /// Common origin: Line Segment & Half Line. Do not keep track of origin
    double commonOrigin(const HalfLine_t& lin, const LineSegment_t& seg, bool backwards=false) const
    { Point_t origin(lin.Start().size()); return commonOrigin(lin, seg, origin, backwards); }
    /// Common origin: Line Segment & Line Segment. Keep track of origin
    double commonOrigin(const HalfLine_t& lin, const LineSegment_t& seg, Point_t& origin, bool backwards=false) const
    { lin.Start().compat(seg.Start()); return _commonOrigin_(lin, seg, origin, backwards); }
    /// Common origin: Line Segment & Half Line. Do not keep track of origin
    double commonOrigin(const LineSegment_t& seg, const HalfLine_t& lin, bool backwards=false) const
    { Point_t origin(lin.Start().size()); return commonOrigin(lin, seg, origin, backwards); }
    /// Common origin: Line Segment & Line Segment. Keep track of origin
    double commonOrigin(const LineSegment_t& seg, const HalfLine_t& lin, Point_t& origin, bool backwards=false) const
    { lin.Start().compat(seg.Start()); return _commonOrigin_(lin, seg, origin, backwards); }
    /// Common origin: Half Line & Half Line. Do not keep track of origin
    double commonOrigin(const HalfLine_t& lin1, const HalfLine_t& lin2, bool backwards=false) const
    { Point_t origin(lin1.Start().size()); return commonOrigin(lin1,lin2, origin, backwards); }
    /// Common origin: Half Line & Half Line. Keep track of origin
    double commonOrigin(const HalfLine_t& lin1, const HalfLine_t& lin2, Point_t& origin, bool backwards=false) const
    { lin1.Start().compat(lin2.Start()); return _commonOrigin_(lin1, lin2, origin, backwards); }
    /// Common origin: Trajectory & Trajectory. Do not keep track of origin
    double commonOrigin(const Trajectory_t& trj1, const Trajectory_t& trj2, bool backwards=false) const
    { Point_t origin(trj1.front().size()); return commonOrigin(trj1,trj2, origin, backwards); }
    /// Common origin: Trajectory & Trajectory. Keep track of origin
    double commonOrigin(const Trajectory_t& trj1, const Trajectory_t& trj2, Point_t& origin, bool backwards=false) const
    { trj1.front().compat(trj2.front()); return _commonOrigin_(trj1, trj2, origin, backwards); }
    /// Common origin: Trajectory & Half Line. Do not keep track of origin
    double commonOrigin(const Trajectory_t& trj, const HalfLine_t& lin, bool backwards=false) const
    { Point_t origin(trj.front().size()); return commonOrigin(trj,lin, origin, backwards); }
    /// Common origin: Trajectory & Half Line. Keep track of origin
    double commonOrigin(const Trajectory_t& trj, const HalfLine_t& lin, Point_t& origin, bool backwards=false) const
    { trj.front().compat(lin.Start()); return _commonOrigin_(trj, lin, origin, backwards); }
    /// Common origin: Trajectory & Half Line. Do not keep track of origin
    double commonOrigin(const HalfLine_t& lin, const Trajectory_t& trj, bool backwards=false) const
    { Point_t origin(trj.front().size()); return commonOrigin(trj,lin, origin, backwards); }
    /// Common origin: Trajectory & Half Line. Keep track of origin
    double commonOrigin(const HalfLine_t& lin, const Trajectory_t& trj, Point_t& origin, bool backwards=false) const
    { trj.front().compat(lin.Start()); return _commonOrigin_(trj, lin, origin, backwards); }
    /// Common origin: Trajectory & Line Segment. Do not keep track of origin
    double commonOrigin(const Trajectory_t& trj, const LineSegment_t& seg, bool backwards=false) const
    { Point_t origin(trj.front().size()); return commonOrigin(trj, seg, origin, backwards); }
    /// Common origin: Trajectory & Line Segment. Keep track of origin
    double commonOrigin(const Trajectory_t& trj, const LineSegment_t& seg, Point_t& origin, bool backwards=false) const
    { trj.front().compat(seg.Start()); return _commonOrigin_(trj, seg, origin, backwards); }
    /// Common origin: Trajectory & Line Segment. Do not keep track of origin
    double commonOrigin(const LineSegment_t& seg, const Trajectory_t& trj, bool backwards=false) const
    { Point_t origin(trj.front().size()); return commonOrigin(trj, seg, origin, backwards); }
    /// Common origin: Trajectory & Line Segment. Keep track of origin
    double commonOrigin(const LineSegment_t& seg, const Trajectory_t& trj, Point_t& origin, bool backwards=false) const
    { trj.front().compat(seg.Start()); return _commonOrigin_(trj, seg, origin, backwards); }

    //************************************************************************
    //BOUNDING SPHERE ALGORITHM: RETURN SMALLEST SPHERE THAT BOUNDS ALL POINTS
    //************************************************************************
    // Bounding Sphere problem given a vector of 3D points
    Sphere_t boundingSphere(const std::vector<Point_t>& pts) const
    { for(auto &p : pts) { pts.front().compat(p); } return _boundingSphere_(pts); }

  protected:

    /// Line & Line distance w/o dimensionality check
    double _SqDist_(const Line_t& l1, const Line_t& l2, Point_t& L1, Point_t& L2) const;

    /// HalfLine & HalfLine distance w/o dimensionality check
    double _SqDist_(const HalfLine_t& l1, const HalfLine_t& l2, Point_t& L1, Point_t& L2) const;

    /// Point & LineSegment distance w/o dimensionality check
    double _SqDist_(const Point_t& pt, const LineSegment_t& line) const
    { return _SqDist_(pt,line.Start(),line.End()); }

    /// Point & LineSegment distance w/o dimensionality check
    double _SqDist_(const Point_t& pt, const Point_t& line_s, const Point_t& line_e) const;

    /// Point & LineSegment distance w/o dimensionality check
    double _SqDist_(const LineSegment_t& line, const Point_t&pt) const
    { return _SqDist_(pt,line); }

    /// HalfLine & LineSegment distance w/o dimensionality check
    double _SqDist_(const HalfLine_t& hline, const LineSegment_t& seg, Point_t& L1, Point_t& L2) const;

    /// LineSegment & LineSegment distance w/o dimensionality check
    double _SqDist_(const LineSegment_t& seg1, const LineSegment_t& seg2, Point_t& c1, Point_t& c2) const;

    // Point & LineSegment closest point w/o dimensionality check
    Point_t _ClosestPt_(const Point_t& pt, const LineSegment_t& line) const;
    // Point & LineSegment closest point w/o dimensionality check
    Point_t _ClosestPt_(const LineSegment_t& line, const Point_t& pt) const
    { return _ClosestPt_(pt,line); }

    /// Point & HalfLine distance w/o dimensionality check
    double _SqDist_(const Point_t& pt, const HalfLine_t& line) const;
    /// Point & HalfLine distance w/o dimensionality check
    double _SqDist_(const HalfLine_t& line, const Point_t& pt) const
    { return _SqDist_(pt,line); }
    // Point & HalfLine closest point w/o dimensionality check
    Point_t _ClosestPt_(const Point_t& pt, const HalfLine_t& line) const;
    // Point & HalfLine closest point w/o dimensionality check
    Point_t _ClosestPt_(const HalfLine_t& line, const Point_t& pt) const
    { return _ClosestPt_(pt,line); }

    // Point & InfLine closest point w/o dimensionality check
    Point_t _ClosestPt_(const Line_t& line, const Point_t& pt) const;
    // Point & InfLine closest point w/o dimensionality check
    Point_t _ClosestPt_(const Point_t& pt, const Line_t& line) const
    { return _ClosestPt_(line,pt); }
    // Point & InfLine  distance w/o dimensionality check
    double _SqDist_(const Line_t& line, const Point_t& pt) const;
    // Point & InfLine  distance w/o dimensionality check
    double _SqDist_(const Point_t& pt, const Line_t& line) const
    { return _SqDist_(line,pt); }

    /// Point & AABox distance w/o dimensionality check
    double _SqDist_(const Point_t& pt, const AABox_t& box) const;
    /// Point & AABox distance w/o dimensionality check
    double _SqDist_(const AABox_t& box, const Point_t& pt) const
    { return _SqDist_(pt,box); }

    /// Point & AABox closest point w/o dimensionality check
    Point_t _ClosestPt_(const Point_t& pt, const AABox_t& box) const;
    /// Point & AABox closest point w/o dimensionality check
    Point_t _ClosestPt_(const AABox_t& box, const Point_t& pt) const
    { return _ClosestPt_(pt,box); }

    /// Common origin: Line & Line. Keep track of origin
    double _commonOrigin_(const Line_t& lin1, const Line_t& lin2, Point_t& origin) const;
    /// Common origin: Half Line & Half Line. Keep track of origin
    double _commonOrigin_(const HalfLine_t& lin1, const HalfLine_t& lin2, Point_t& origin, bool backwards) const;
    /// Common origin: Line Segment & Half Line. Keep track of origin
    double _commonOrigin_(const HalfLine_t& lin, const LineSegment_t& seg, Point_t& origin, bool backwards) const;
    /// Common origin: Line Segment & Line Segment. Keep track of origin
    double _commonOrigin_(const LineSegment_t& seg1, const LineSegment_t& seg2, Point_t& origin, bool backwards) const;
    /// Common origin: Trajectory & Trajectory. Keep track of origin
    double _commonOrigin_(const Trajectory_t& trj1, const Trajectory_t& trj2, Point_t& origin, bool backwards) const;
    /// Common origin: Trajectory & Line Segment. Keep track of origin
    double _commonOrigin_(const Trajectory_t& trj, const LineSegment_t& seg, Point_t& origin, bool backwards) const;
    /// Common origin: Trajectory & Half Line. Keep track of origin
    double _commonOrigin_(const Trajectory_t& trj, const HalfLine_t& lin, Point_t& origin, bool backwards) const;


    // Bounding Sphere given a vector of points
    Sphere_t _boundingSphere_(const std::vector<Point_t>& pts) const;
    Sphere_t _RemainingPoints_(std::vector<Point_t>& remaining, const Sphere_t& thisSphere)const;
    Sphere_t _WelzlSphere_(const std::vector<Point_t>& pts, int numPts, std::vector<Point_t> sosPts) const;

    /// Clamp function: checks if value out of bounds
    double _Clamp_(const double n, const double min, const double max) const;

    /// Swap two points if min & max are inverted
    inline void _Swap_(double& tmin, double& tmax) const
    { if(tmin > tmax) std::swap(tmin,tmax); }

  };
}

#endif
/** @} */ // end of doxygen group 

