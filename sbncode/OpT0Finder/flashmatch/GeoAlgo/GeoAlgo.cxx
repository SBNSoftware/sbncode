#ifndef BASICTOOL_GEOALGO_CXX
#define BASICTOOL_GEOALGO_CXX

#include "GeoAlgo.h"
#include <cassert>

namespace geoalgo {

  // Ref. RTCD 5.3.2 p. 177
  // Intersection of a HalfLine w/ AABox
  std::vector<Point_t> GeoAlgo::Intersection(const AABox_t& box,
					     const HalfLine_t& line,
					     bool back) const
  {
    // Result container
    std::vector<Point_t> result;
    Point_t xs1(3); // Only 2 points max possible
    Point_t xs2(3); // Create in advance for early termination checks
    // One-time only initialization for unit vectors
    static std::vector<Point_t> min_plane_n;
    static std::vector<Point_t> max_plane_n;
    if(!min_plane_n.size()) {
      min_plane_n.reserve(3);
      min_plane_n.push_back(Vector_t(-1,0,0));
      min_plane_n.push_back(Vector_t(0,-1,0));
      min_plane_n.push_back(Vector_t(0,0,-1));
      max_plane_n.reserve(3);
      max_plane_n.push_back(Vector_t(1,0,0));
      max_plane_n.push_back(Vector_t(0,1,0));
      max_plane_n.push_back(Vector_t(0,0,1));
    }
    // min/max points of the AABox
    auto const& min_pt = box.Min();
    auto const& max_pt = box.Max();
    // start/dir of the line
    auto const& start = line.Start();
    auto        dir   = line.Dir();
    if(back) dir *= -1;
    // Inspect the case of parallel line
    for(size_t i=0; i<min_pt.size(); ++i) {
      if ( dir[i] == 0 &&
	   (start[i] <= min_pt[i] ||  max_pt[i] <= start[i]) )
	return result;
    }
    // Look for xs w/ 3 planes
    for(size_t i=0; i<3; ++i) {
      auto const& normal = min_plane_n[i];
      double s = (-1.) * ( normal * (start - min_pt) ) / (normal * dir);
      if(s<0) continue;
      auto xs = start + dir * s;
      // Check if the found point is within the surface area of other 2 axis
      bool on_surface=true;
      for(size_t sur_axis=0; sur_axis<3; ++sur_axis) {
	if(sur_axis==i) continue;
	if(xs[sur_axis] < min_pt[sur_axis] || max_pt[sur_axis] < xs[sur_axis]) {
	  on_surface=false;
	  break;
	}
      }
      if(on_surface && xs != xs1) {
	// Directly assign to xs1 instead of making a copy
	if(!(xs1.IsValid())) for(size_t j=0; j<3; ++j) xs1[j]=xs[j];
	else {
	  // If xs2 is filled, no more point to search. Exit.
	  for(size_t j=0; j<3; ++j) xs2[j]=xs[j];
	  break;
	}
      }
    }
    // If xs2 is filled, simply return the result. Order the output via distance
    if(xs2.IsValid()) {
      result.reserve(2);
      if(xs1._SqDist_(start) > xs2._SqDist_(start)) std::swap(xs1,xs2);
      result.push_back(xs1);
      result.push_back(xs2);
      return result;
    }
    // Look for xs w/ 3 planes
    for(size_t i=0; i<3; ++i) {
      auto const& normal = max_plane_n[i];
      double s = (-1.) * ( normal * (start - max_pt) ) / (normal * dir);
      if(s<0) continue;
      auto xs = start + dir * s;
      bool on_surface=true;
      for(size_t sur_axis=0; sur_axis<3; ++sur_axis) {
	if(sur_axis==i) continue;
	if(xs[sur_axis] < min_pt[sur_axis] || max_pt[sur_axis] < xs[sur_axis]) {
	  on_surface=false;
	  break;
	}
      }
      if(on_surface && xs != xs1) {
	if(!(xs1.IsValid())) for(size_t j=0; j<3; ++j) xs1[j]=xs[j];
	else {
	  for(size_t j=0; j<3; ++j) xs2[j]=xs[j];
	  break;
	}
      }
    }
    if(!xs1.IsValid()) return result;
    if(xs2.IsValid()) {
      result.reserve(2);
      if(xs1._SqDist_(start) > xs2._SqDist_(start)) std::swap(xs1,xs2);
      result.push_back(xs1);
      result.push_back(xs2);
      return result;
    }
    result.push_back(xs1);
    return result;
  }

  // AABox_t & LineSegment_t intersection search. Make a use of AABox_t & HalfLine_t function
  std::vector<Point_t> GeoAlgo::Intersection(const AABox_t& box,
					     const LineSegment_t& line) const
  {
    auto const& st = line.Start();
    auto const& ed = line.End();
    // Create a static HalfLine_t for this algorithm
    static HalfLine_t hline(st,ed-st);
    hline.Start(st[0],st[1],st[2]);
    hline.Dir(ed[0]-st[0],ed[1]-st[1],ed[2]-st[2]);
    auto hline_result = Intersection(box,hline);

    if(!hline_result.size()) return hline_result;

    // For non-empty result, only keep ones that is within the line length
    std::vector<Point_t> result;
    auto length = st._SqDist_(ed);
    for(auto const& pt : hline_result)
      if(st._SqDist_(pt) < length) result.push_back(pt);
    return result;
  }

  // AABox_t & Trajectory_t intersection search. Make a use of AABox_t & HalfLine_t function
  std::vector<Point_t> GeoAlgo::Intersection(const AABox_t& box,
					     const Trajectory_t& trj) const
  {

    std::vector<Point_t> result;
    if(trj.size() < 2) return result; // If only 1 point, return
    // Check compat
    trj.compat(box.Min());
    // Call the onetime-only HalfLine constructor
    static HalfLine_t hline(trj[0],trj[1]);
    for(size_t i=0; i<trj.size()-1; ++i) {

      auto const& st = trj[i];
      auto const& ed = trj[i+1];
      hline.Start(st[0],st[1],st[2]);
      hline.Dir(ed[0]-st[0],ed[1]-st[1],ed[2]-st[2]);

      auto hline_result = Intersection(box,hline);

      if(!hline_result.size()) continue;

      // Check if the length makes sense
      auto length = st._SqDist_(ed);
      for(auto const& pt : hline_result)
	if(st._SqDist_(pt) < length) result.push_back(pt);

    }
    return result;
  }

  // LineSegment sub-segment of HalfLine inside an AABox w/o checks
  LineSegment_t GeoAlgo::BoxOverlap(const AABox_t& box, const HalfLine_t& line) const
  {
    // First find interection point of half-line and box
    auto xs_v = Intersection(box, line);
    if(xs_v.size()==2) return LineSegment_t(xs_v[0],xs_v[1]);

    // Build a new LineSegment
    if(!(xs_v.size())) return LineSegment_t();

    // Only other possiblity is # = 1
    return LineSegment_t(line.Start(),xs_v[0]);
  }

  /// Get Trajectory inside box given some input trajectory -> now assumes trajectory cannot exit and re-enter box
  Trajectory_t GeoAlgo::BoxOverlap(const AABox_t& box, const Trajectory_t& trj) const
  {

    // if empty, return as is
    if(trj.empty()) return trj;
    // if first & last points inside, then return full trajectory
    if ( box.Contain(trj[0]) && box.Contain(trj.back()) )
      return trj;

    Trajectory_t result;
    result.reserve(trj.size());
    // This assumes points are ordered along the trajectory.
    for(size_t i=0; (i+1)<trj.size(); ++i) {
      auto const& pt0 = trj[i];
      auto const& pt1 = trj[i+1];
      if(box.Contain(pt0)) result.push_back(pt0);
      // Are we stepping through the boundary
      // need to account for both directions of trajectory
      if ((box.Contain(pt0) && !box.Contain(pt1)) || (!box.Contain(pt0) && box.Contain(pt1))) {
            // to make sure we find the intersection even if one of the points is exactly
            // at the boundary, displace a little bit both points along the direction
            geoalgo::Vector ptA(pt0[0], pt0[1], pt0[2]);
            geoalgo::Vector ptB(pt1[0], pt1[1], pt1[2]);
            geoalgo::Vector unit = (ptB - ptA).Dir();
            geoalgo::Vector small_unit(unit[0]*0.0001, unit[1]*0.0001, unit[2]*0.0001);
            ptA = ptA - small_unit;
            ptB = ptB + small_unit;
            // std::cout << ptA << ptB << std::endl;
            geoalgo::LineSegment seg(ptA[0], ptA[1], ptA[2], ptB[0], ptB[1], ptB[2]);

            // now find intersecting point
        	auto crossing_pts = this->Intersection(box,seg);

            for (auto pt : crossing_pts) {
                // std::cout << pt << " " << box.Contain(pt) << std::endl;
                if (box.Contain(pt)) result.push_back(pt);
            }
            // std::cout << crossing_pts.size() << std::endl;
            // There can be 2 intersection points in very rare cases....
            // Pt (-368.490041,129.954320,-587.387966) Pt (-366.214146,134.960090,-588.246902)
            // this segment intersects both x and y planes
        	// assert(crossing_pts.size()==1);
        	// result.push_back(crossing_pts[0]);
    }
    }

    if(box.Contain(trj.back()))
      result.push_back(trj.back());
    // else if(trj.size()>1) {
    //   auto const& pt0 = trj[trj.size()-2];
    //   auto const& pt1 = trj.back();
    //   if(box.Contain(pt0)) {
	// geoalgo::LineSegment seg(pt0,pt1);
	// auto crossing_pts = this->Intersection(box,seg);
	// assert(crossing_pts.size()==1);
	// result.push_back(crossing_pts[0]);
    //   }
    // }

    return result;
  }

  // Ref. RTCD 5.1.8 p. 146
  // Distance between two infinite lines
  double GeoAlgo::_SqDist_(const Line_t& l1, const Line_t &l2, Point_t& L1, Point_t& L2) const
  {

    // closest approach when segment connecting the two lines
    // is perpendicular to both lines

    // L1 = P1 + s(Q1-P1)
    // L1 = P2 + t(Q2-P2)
    // L1(s) and L2(t) are the closest approach points
    // d1 = Q1-P1
    // d2 = Q2-P2
    // v(s,t) = L1(s) - L2(t)
    // require d1*v == 0 && d2*v == 0
    Vector_t d1 = l1.Pt2()-l1.Pt1();
    Vector_t d2 = l2.Pt2()-l2.Pt1();
    Vector_t r = l1.Pt1()-l2.Pt1();

    double a = d1*d1;
    double b = d1*d2;
    double c = d1*r;
    double e = d2*d2;
    double f = d2*r;

    double d = a*e-b*b;

    // if d==0 the lines are parallel
    // return Pt1 (doesn't matter) for line 1
    // distance is distance between Pt1 of 1 line & second line
    if (d == 0){
      L1 = l1.Pt1();
      L2 = _ClosestPt_(l1.Pt1(),l2);
      return L1._SqDist_(L2);
    }

    double s = (b*f-c*e)/d;
    double t = (a*f-b*c)/d;

    // s & t represent the paramteric points on the lines
    // for the closest approach point
    // now find the Point_t object at those locations

    L1 = l1.Pt1() + ( l1.Pt2()-l1.Pt1() )*s;
    L2 = l2.Pt1() + ( l2.Pt2()-l2.Pt1() )*t;

    // find distance between these points
    double dist = L1._SqDist_(L2);

    return dist;
  }


  // Distance between two half-infinite lines
  // use same function as for infinite lines
  // but expect return points L1 and L2 to be
  // "after" start point. Otherwise need
  // to re-calculate using start point
  // for one or both of the half-lines
  double GeoAlgo::_SqDist_(const HalfLine_t& l1, const HalfLine_t &l2, Point_t& L1, Point_t& L2) const
  {

    //Same as for _SqDist_ with infinite line but check whether s & t go out of bounds (i.e. negative)

    Vector_t d1 = l1.Dir();
    Vector_t d2 = l2.Dir();
    Vector_t r = l1.Start()-l2.Start();

    double a = d1*d1;
    double b = d1*d2;
    double c = d1*r;
    double e = d2*d2;
    double f = d2*r;

    double d = a*e-b*b;

    // Need to make sure d != 0
    if ( d==0 ){
      // lines are parallel
      // check closest distance from one start point
      // to other line. Order indifferent
      // Choose l1 to have cloest point be Start point
      L1 = l1.Start();
      L2 = _ClosestPt_(L1,l2);
      return L1._SqDist_(L2);
    }

    double s = (b*f-c*e)/d;
    double t = (a*f-b*c)/d;

    // if s or t < 0, out of bounds for half-line
    if (s < 0) s = 0;
    if (t < 0) t = 0;

    // s & t represent the paramteric points on the lines
    // for the closest approach point
    // now find the Point_t object at those locations

    L1 = l1.Start() + l1.Dir()*s;
    L2 = l2.Start() + l2.Dir()*t;

    // find distance between these points
    double dist = L1._SqDist_(L2);

    return dist;
  }



  // Distance between two half-infinite lines
  // use same function as for infinite lines
  // but expect return points L1 and L2 to be
  // "after" start point. Otherwise need
  // to re-calculate using start point
  // for one or both of the half-lines
  double GeoAlgo::_SqDist_(const HalfLine_t& hline, const LineSegment_t &seg, Point_t& L1, Point_t& L2) const
  {

    //Same as for _SqDist_ with infinite line but check whether s & t go out of bounds (i.e. negative)

    auto const  d1 = hline.Dir();
    auto const  d2 = seg.End()-seg.Start();
    auto const   r = hline.Start()-seg.Start();

    double a = d1*d1;
    double b = d1*d2;
    double c = d1*r;
    double e = d2*d2;
    double f = d2*r;

    double d = a*e-b*b;

    // if parallel then d == 0
    if (d == 0){
      // distance is smallest quantity between:
      // - distance from segment start to line
      // - distance from segment end to line
      double sDist = _SqDist_(seg.Start(),hline);
      double eDist = _SqDist_(seg.End(),hline);
      if ( sDist <= eDist ) {
	// get closest point on line
	L1 = _ClosestPt_(seg.Start(),hline);
	L2 = seg.Start();
	return sDist;
      }
      else{
	L1 = _ClosestPt_(seg.End(),hline);
	L2 = seg.End();
	return eDist;
      }
    }// if parallel

    double s = (b*f-c*e)/d;


    // now check if parameters are out of bounds
    if ( s < 0 ){
      s = 0; // closest point on half-line is start
      // re-evaluate closest point on segment using line start point
      L1 = hline.Start();
      L2 = _ClosestPt_(L1,seg);
      return L1._SqDist_(L2);
    }

    // if closest point is not beyond half-line
    // it could be due to an intersection between half-line
    // and segment projection.
    // check value of t
    double t = (a*f-b*c)/d;
    // if t > 0 && < 1 then the two lines intersect. We are good!
    if ( (t < 1) && (t > 0) ){
      L1 = hline.Start() + hline.Dir()*s;
      L2 = seg.Start() + (seg.End()-seg.Start())*t;
      return L1._SqDist_(L2);
    }
    // if out of bounds clamp
    // then re-evaluate closest point on line
    t = _Clamp_(t,0,1);
    L2 = seg.Start() + (seg.End()-seg.Start())*t;
    L1 = _ClosestPt_(L2,hline);
    return L1._SqDist_(L2);
  }

  // Ref. RTCD Ch 5.1 p. 130
  double GeoAlgo::_SqDist_(const Point_t& pt, const Point_t& line_s, const Point_t& line_e) const
  {
    auto const  ab = line_e - line_s;
    auto const  ac = pt - line_s;
    auto const  bc = pt - line_e;
    auto e = ac * ab;
    if( e <= 0. ) return ac.SqLength();
    auto f = ab.SqLength();
    if( e >= f ) return bc.SqLength();
    return (ac.SqLength() - e * e / f);
  }

  // Ref. RTCD Ch 5.1 p. 128-129
  Point_t GeoAlgo::_ClosestPt_(const Point_t& pt, const LineSegment_t& line) const
  {

    auto const& ab = line.Dir();
    // Project pt on line (ab), but deferring divide by ab * ab
    auto t = ((pt - line.Start()) * ab);
    // pt projects outside line, on the start side; clamp to start
    if( t <= 0. ) return line.Start();
    else {
      auto denom = ab.SqLength();
      // pt projects outside line, on the end side; clamp to end
      if( t>= denom ) return line.End();
      // pt projects inside the line. must deferred divide now
      else return (line.Start() + ab * (t/denom));
    }
  }

  // Ref. RTCD Ch 5.1 p. 130
  double GeoAlgo::_SqDist_(const Point_t& pt, const HalfLine_t& line) const
  {
    auto const& ab = line.Dir();
    auto const  ac = pt - line.Start();
    auto const  bc = ac - ab;

    auto e = ac * ab;
    if( e <= 0. ) return (ac * ac);
    auto f = ab.SqLength();
    return (ac.SqLength() - e * e / f);
  }


  // Ref. RTCD Ch 5.1 p. 128-129
  Point_t GeoAlgo::_ClosestPt_(const Point_t& pt, const HalfLine_t& line) const
  {
    auto const& ab = line.Dir();
    auto t = (pt - line.Start()) * ab;
    if( t <= 0. ) return line.Start();
    else {
      auto denom = ab.Length();
      return (line.Start() + ab * (t/denom));
    }
  }

  // Point & Infinite Line min Distance
  double GeoAlgo::_SqDist_(const Line_t& line, const Point_t& pt) const
  {
    auto const  ab = line.Pt2()-line.Pt1();
    auto const  ac = pt - line.Pt1();
    auto const  bc = ac - ab;

    auto e = ac * ab;
    auto f = ab.SqLength();
    return (ac.SqLength() - e * e / f);
  }

  // Point & Infinite Line Closest Point
  Point_t GeoAlgo::_ClosestPt_(const Line_t& line, const Point_t& pt) const
  {
    auto const& ab = line.Pt2()-line.Pt1();
    auto t = (pt - line.Pt1()) * ab;
    auto denom = ab.SqLength();
    return (line.Pt1() + ab * (t/denom));
  }


  // Ref. RTCD Ch 5.1 p. 131-132 ... modified to consider distance to the box's wall
  double GeoAlgo::_SqDist_(const Point_t& pt, const AABox_t& box) const
  {
    double dist = kINVALID_DOUBLE;

    // If a point is inside the box, simply compute the smallest perpendicular distance
    if(box.Contain(pt)) {

      auto const& pt_min = box.Min();
      auto const& pt_max = box.Max();
      // (1) Compute the distance to the YZ wall
      double dist_to_yz = pt[0] - pt_min[0];
      if( dist_to_yz > (pt_max[0] - pt[0]) ) dist_to_yz = pt_max[0] - pt[0];

      // (2) Compute the distance to the XZ wall
      double dist_to_zx = pt[1] - pt_min[1];
      if( dist_to_zx > (pt_max[1] - pt[1]) ) dist_to_zx = pt_max[1] - pt[1];

      // (3) Compute the distance to the XY wall
      double dist_to_xy = pt[2] - pt_min[2];
      if( dist_to_xy > (pt_max[2] - pt[2]) ) dist_to_xy = pt_max[2] - pt[2];

      // (4) Compute the minimum of (3), (4), and (5)
      dist = ( dist_to_yz < dist_to_zx ? dist_to_yz : dist_to_zx );
      dist = ( dist < dist_to_xy ? dist : dist_to_xy );
      dist *= dist;

    }

    else{
      // This refers to Ref. RTCD 5.1.3.1
      // re-set distance
      dist = 0;
      for(size_t i=0; i<pt.size(); ++i) {

	auto const& v_pt  = pt[i];
	auto const& v_max = box.Max()[i];
	auto const& v_min = box.Min()[i];

	if(v_pt < v_min) dist += (v_min - v_pt) * (v_min - v_pt);
	if(v_pt > v_max) dist += (v_pt - v_max) * (v_pt - v_max);
      }
    }
    return dist;
  }

  // Ref. RTCD Ch 5.1 p. 130-131 ... modified to consider a point on the surface
  Point_t GeoAlgo::_ClosestPt_(const Point_t& pt, const AABox_t& box) const
  {
    // Result
    auto res = pt;
    // For each coordinate axis, if the point coordinate value is outside box,
    // clamp it to the box, else keep it as is
    for(size_t i=0; i<pt.size(); ++i) {
      auto const& v_pt  = pt[i];
      auto const& v_min = box.Min()[i];
      auto const& v_max = box.Max()[i];
      res[i] = v_pt;
      if( v_pt < v_min ) res[i] = v_min;
      if( v_pt > v_max ) res[i] = v_max;
    }
    return res;
  }



  // Distance between a Trajectory_t and a Point_t
  // Loop over segments that make up the trajectory and keep track
  // of shortest distance between any of them and the point
  // The smallest such distance is the return
  double GeoAlgo::SqDist(const Point_t& pt, const Trajectory_t& trj) const
  {

    // Make sure trajectory object is properly defined
    if (!trj.size())
      throw GeoAlgoException("Trajectory object not properly set...");

    // Check dimensionality compatibility between point and trajectory
    trj.compat(pt);

    // Now keep track of smallest distance and loop over traj segments
    double distMin = kINVALID_DOUBLE;
    for (size_t l=0; l < trj.size()-1; l++){
      double distTmp = _SqDist_(pt,trj[l],trj[l+1]);
      if (distTmp < distMin) { distMin = distTmp; }
    }

    return distMin;
  }



  // Distance between vector of Trajectories and a Point
  // Loop over Trajectories and find the closest one
  // then keep track of that closest one
  double GeoAlgo::SqDist(const Point_t& pt, const std::vector<Trajectory_t> &trj, int &trackIdx) const
  {

    // holder for shortest distance
    double minDist = kINVALID_DOUBLE;

    // loop over trajectories
    for (size_t t=0; t < trj.size(); t++){

      // trajectory & point dimensionality will be checked in next function
      // now calculate distance w.r.t. this track
      double tmpDist = SqDist(pt, trj[t]);
      if (tmpDist < minDist){
	minDist = tmpDist;
	trackIdx = t;
      }
    }// for all tracks

    return minDist;
  }



  // Closest point between a Trajectory and a Point
  // Loop over segments that make up the trajectory and keep track
  // of shortest distance between any of them and the point
  // For the shortest distance find the point at which it is found
  Point_t GeoAlgo::ClosestPt(const Point_t& pt, const Trajectory_t& trj, int &idx) const
  {

    // Make sure trajectory object is properly defined
    if (!trj.size())
      throw GeoAlgoException("Trajectory object not properly set...");

    // Check dimensionality compatibility between point and trajectory
    trj.compat(pt);

    // Now keep track of smallest distance and loop over traj segments
    double distMin = kINVALID_DOUBLE;
    // For that smallest distance, keep track of the segment for which it was found
    for (size_t l=0; l < trj.size()-1; l++){
      double distTmp = _SqDist_(pt,trj[l],trj[l+1]);
      if (distTmp < distMin) { distMin = distTmp; idx = l; }
    }

    // Now that we have the segment for the closest approach
    // Use it to find the closest point on that segment
    LineSegment_t segMin(trj[idx], trj[idx+1]);
    return _ClosestPt_(pt,segMin);
  }


  // Closest point between a vector of trajectories and a point
  // Loop over segments that make up the trajectory and keep track
  // of shortest distance between any of them and the point
  // For the shortest distance find the point at which it is found
  Point_t GeoAlgo::ClosestPt(const Point_t& pt, const std::vector<Trajectory_t> &trj, int& TrackIdx) const
  {

    // since finding the smallest distance is faster than finding the closest point
    // first loop through tracks, and find the one that is closest to the point
    // then finally find the closest point on that track

    for (size_t t=0; t < trj.size(); t++){

      // holder for smallest distance
      double minDist = kINVALID_DOUBLE;

      // track & point dimensionality will be checked per-track by next function
      // now calculate distance w.r.t. this track
      double tmpDist = SqDist(pt, trj[t]);
      if (tmpDist < minDist){
	minDist = tmpDist;
	TrackIdx = t;
      }

    }// for all tracks

    // now we know which track is closest -> find the closest point to that track
    return ClosestPt(pt, trj[TrackIdx]);
  }


  // Closest Approach between a segment and a Trajectory
  // loop over segments in trajectory and find the one that
  // is closest. Then find distance
  double GeoAlgo::SqDist(const LineSegment_t& seg, const Trajectory_t& trj, Point_t& c1, Point_t& c2) const
  {

    // Make sure trajectory object is properly defined
    if (!trj.size())
      throw GeoAlgoException("Trajectory object not properly set...");

    // Check dimensionality compatibility between point and trajectory
    trj.compat(seg.Start());

    // keep track of c1 & c2
    Point_t c1min;
    Point_t c2min;
    // Now keep track of smallest distance and loop over traj segments
    double distMin = kMAX_DOUBLE;

    for (size_t l=0; l < trj.size()-1; l++){
      LineSegment_t segTmp(trj[l], trj[l+1]);
      double distTmp = _SqDist_(segTmp, seg, c1min, c2min);
      if ( distTmp < distMin ){
	c1 = c1min;
	c2 = c2min;
	distMin = distTmp;
      }
    }//for all segments in the track

    return distMin;
  }


  // Closest Approach between a Trajectory and a Trajectory
  // loop over segments in trajectory1 and those in
  // trahectory2 and find the best point
  double GeoAlgo::SqDist(const Trajectory_t& trj1, const Trajectory_t& trj2, Point_t& c1, Point_t& c2) const
  {

    // Make sure trajectory object is properly defined
    if ( !trj1.size() or !trj2.size() )
      throw GeoAlgoException("Trajectory object not properly set...");

    // Check dimensionality compatibility between point and trajectory
    trj1.compat(trj2[0]);

    // keep track of c1 & c2
    Point_t c1min;
    Point_t c2min;
    // Now keep track of smallest distance and loop over traj segments
    double distMin = kMAX_DOUBLE;

    for (size_t l1=0; l1 < trj1.size()-1; l1++){
      for (size_t l2=0; l2 < trj2.size()-1; l2++){
	LineSegment_t segTmp1(trj1[l1], trj1[l1+1]);
	LineSegment_t segTmp2(trj2[l2], trj2[l2+1]);
	double distTmp = _SqDist_(segTmp1, segTmp2, c1min, c2min);
	if ( distTmp < distMin ){
	  c1 = c1min;
	  c2 = c2min;
	  distMin = distTmp;
	}
      }// for segments in trajectory 2
    }//for all segments in trajectory 1

    return distMin;
  }

  // Closest Approach between a HalfLine and a Trajectory
  // loop over segments in trajectory and find the one that
  // is closest. Then find distance
  double GeoAlgo::SqDist(const HalfLine_t& hline, const Trajectory_t& trj, Point_t& c1, Point_t& c2) const
  {

    // Make sure trajectory object is properly defined
    if (!trj.size())
      throw GeoAlgoException("Trajectory object not properly set...");

    // Check dimensionality compatibility between point and trajectory
    trj.compat(hline.Start());

    // keep track of c1 & c2
    Point_t c1min;
    Point_t c2min;
    // Now keep track of smallest distance and loop over traj segments
    double distMin = kMAX_DOUBLE;

    for (size_t l=0; l < trj.size()-1; l++){
      LineSegment_t segTmp(trj[l], trj[l+1]);
      double distTmp = _SqDist_(hline, segTmp, c1min, c2min);
      if ( distTmp < distMin ){
	c1 = c1min;
	c2 = c2min;
	distMin = distTmp;
      }
    }//for all segments in the track

    return distMin;
  }


  // Closest Approach between a Segment and a vector of tracks
  // keep track of points of closest approach both on nearest
  // track as well as on segment
  // keep track of which track has the point of closest approcah
  double GeoAlgo::SqDist(const LineSegment& seg, const std::vector<Trajectory_t> &trj, Point_t& c1, Point_t& c2, int& trackIdx) const
  {

    // holders to keep track of track with shortest distance
    double minDist = kMAX_DOUBLE;
    // holders for points of closest approach
    Point_t c1min;
    Point_t c2min;

    for (size_t t=0; t < trj.size(); t++){

      //need to loop over all tracks and find the one which is closest

      // dimensionality checks will be done in next function, per track.

      // now calculate closest approach
      double tmpDist = SqDist(seg, trj[t], c1min, c2min);

      // is this the best yet?
      if (tmpDist < minDist){
	minDist = tmpDist;
	c1 = c1min;
	c2 = c2min;
	trackIdx = t;
      }

    }// for all tracks in vector

    return minDist;
  }

  // Ref. RTCD Sec. 5.1.9 - pg. 148-150
  double GeoAlgo::_SqDist_(const LineSegment_t& seg1, const LineSegment_t& seg2,
				Point_t& c1, Point_t& c2) const
  {

    double t1, t2;

    auto const& s1 = seg1.Start();
    auto const& s2 = seg2.Start();
    auto const& e1 = seg1.End();
    auto const& e2 = seg2.End();

    auto d1 = e1 - s1;
    auto d2 = e2 - s2;
    auto r  = s1 - s2;

    double a = d1.SqLength();
    double e = d2.SqLength();
    double f    = d2 * r;

    // check if segment is too short
    if ( (a <= 0) && (e <= 0) ){
      //both segments are too short
      t1 = t2 = 0.;
      c1 = s1;
      c2 = s2;
      Vector_t distVector = c2 - c1;
      return distVector.SqLength();
    }
    if (a <= 0){
      //first segment degenerates into a point
      t1 = 0.;
      t2 = f/e;
      t2 = _Clamp_(t2,0.,1.);
    }
    else{
      double c = d1 * r;
      if (e <= 0){
	//second segment degenerates into a point
	t2 = 0.;
	t1 = _Clamp_(-c/a,0.,1.);
      }
      else{
	// the general case...no degeneracies
	double b = d1 * d2;
	double denom = (a*e)-(b*b);

	if (denom != 0.)
	  t1 = _Clamp_((b*f-c*e)/denom, 0., 1.);
	else
	  t1 = 0.;

	t2 = (b*t1+f)/e;

	if (t2 < 0.){
	  t2 = 0.;
	  t1 = _Clamp_(-c/a, 0., 1.);
	}
	else if (t2 > 1.){
	  t2 = 1.;
	  t1 = _Clamp_((b-c)/a, 0., 1.);
	}

      }
    }

    c1 = s1 + d1 * t1;
    c2 = s2 + d2 * t2;

    Vector_t distVector = c2 - c1;
    return distVector.SqLength();

  }


  // Clamp function:
  // if 1st argument out of bounds w.r.t. min & max
  // return the boundary point
  double GeoAlgo::_Clamp_(const double n, const double min, const double max) const
  {
    if (n < min) { return min; }
    if (n > max) { return max; }
    return n;
  }


  /// Common origin: Half Line & Half Line. Keep track of origin
  double GeoAlgo::_commonOrigin_(const Line_t& lin1, const Line_t& lin2, Point_t& origin) const
  {

    // Function Description:
    // Given two HalfLine objects, project them backwards
    // and find the point of closest approach on both
    // lines.
    // The half-point between these two is considered the
    // candidate origin or vertex of the two lines
    // Then make segments uniting this vertex and
    // the start point of both half lines.
    // Take the dot product of the segment uniting
    // the vertex with the start of lin1, and the direction
    // of lin1. Similarly for lin2.
    // These dot-products will be close to 1 if the lines,
    // traced backwards, indeed point to the reconstructed
    // vertex.
    // return the sum of these dot products, which
    // is bound between -2 and +2.
    // other values of this return (1, 0, -1, -2)
    // will give insight on other possible topologies.

    // get directions of two lines
    Vector_t dir1(lin1.Pt2()-lin1.Pt1());
    Vector_t dir2(lin2.Pt2()-lin2.Pt1());
    dir1.Normalize();
    dir2.Normalize();

    // Closest approach points on the two lines
    Point_t pt1(lin1.Pt1().size());
    Point_t pt2(lin2.Pt1().size());

    //double IP = _SqDist_(lin1, lin2, pt1, pt2);
    origin = (pt1+pt2)/2.;

    // If origin coincides with lin1 start
    // -> vec1 should be in same direction of lin1
    Vector_t vec1(dir1);
    if (lin1.Pt1() != origin)
      vec1 = lin1.Pt1()-origin;
    vec1.Normalize();
    // similarly for vec1
    Vector_t vec2(dir2);
    if (lin2.Pt1() != origin)
      vec2 = lin2.Pt1()-origin;
    vec2.Normalize();

    return vec1.Dot(dir1) + vec2.Dot(dir2);
  }


  /// Common origin: Half Line & Half Line. Keep track of origin
  double GeoAlgo::_commonOrigin_(const HalfLine_t& lin1, const HalfLine_t& lin2, Point_t& origin, bool backwards) const
  {

    //If backwards is false, call infinite line function, otherwise proceed
    if (!backwards){
      Line_t l1(lin1.Start(),lin1.Start()+lin1.Dir());
      Line_t l2(lin2.Start(),lin2.Start()+lin2.Dir());
      return _commonOrigin_(l1,l2,origin);
    }

    // Function Description:
    // Given two HalfLine objects, project them backwards
    // and find the point of closest approach on both
    // lines.
    // The half-point between these two is considered the
    // candidate origin or vertex of the two lines
    // Then make segments uniting this vertex and
    // the start point of both half lines.
    // Take the dot product of the segment uniting
    // the vertex with the start of lin1, and the direction
    // of lin1. Similarly for lin2.
    // These dot-products will be close to 1 if the lines,
    // traced backwards, indeed point to the reconstructed
    // vertex.
    // return the sum of these dot products, which
    // is bound between -2 and +2.
    // other values of this return (1, 0, -1, -2)
    // will give insight on other possible topologies.

    // Flip the HalfLines: want to project backwards
    HalfLine_t lin1Back(lin1.Start(), lin1.Dir()*(-1));
    HalfLine_t lin2Back(lin2.Start(), lin2.Dir()*(-1));
    // Closest approach points on the two lines
    Point_t pt1(lin1.Start().size());
    Point_t pt2(lin2.Start().size());

    _SqDist_(lin1Back, lin2Back, pt1, pt2);
    origin = (pt1+pt2)/2.;

    // If origin coincides with lin1 start
    // -> vec1 should be in same direction of lin1
    Vector_t vec1(lin1.Dir());
    if (lin1.Start() != origin)
      vec1 = lin1.Start()-origin;
    vec1.Normalize();
    // similarly for vec1
    Vector_t vec2(lin2.Dir());
    if (lin2.Start() != origin)
      vec2 = lin2.Start()-origin;
    vec2.Normalize();

    return vec1.Dot(lin1.Dir()) + vec2.Dot(lin2.Dir());
    //    std::cout << "dot is: "  << dot << std::endl;
    //    if ( !((dot <= 2) && (dot >= -2)) )
    //      throw GeoAlgoException("commonOrigin failed. Sum of two dot-products must be bound by [-2,2]");

    //    return dot;
  }
  /// Common origin: Line Segment & Half Line. Keep track of origin
  double GeoAlgo::_commonOrigin_(const HalfLine_t& lin, const LineSegment_t& seg, Point_t& origin, bool backwards) const
  {
    // Make a Half-line out of the line-segment
    // we want to project backwards to a common origin
    // not limit ourselves to an origin that must be on the segment
    HalfLine_t lin2(seg.Start(), seg.Dir());
    return _commonOrigin_(lin, lin2, origin, backwards);
  }
  /// Common origin: Line Segment & Line Segment. Keep track of origin
  double GeoAlgo::_commonOrigin_(const LineSegment_t& seg1, const LineSegment_t& seg2, Point_t& origin, bool backwards) const
  {
    // Make a Half-line out of the line-segments
    // we want to project backwards to a common origin
    // not limit ourselves to an origin that must be on the segment
    HalfLine_t lin1(seg1.Start(), seg1.Dir());
    HalfLine_t lin2(seg2.Start(), seg2.Dir());
    return _commonOrigin_(lin1, lin2, origin, backwards);
  }

  /// Common origin: Trajectory & Trajectory/ Keep track of origin
  double GeoAlgo::_commonOrigin_(const Trajectory_t& trj1, const Trajectory_t& trj2, Point_t& origin, bool backwards) const
  {
    // Turn the trajectory into half-line that connect start -> end
    HalfLine_t lin1(trj1.front(),trj1.back()-trj1.front());
    // Turn the segment into half-line
    HalfLine_t lin2(trj2.front(),trj2.back()-trj2.front());
    return _commonOrigin_(lin1, lin2, origin, backwards);
  }

  /// Common origin: Trajectory & Line Segment. Keep track of origin
  double GeoAlgo::_commonOrigin_(const Trajectory_t& trj, const LineSegment_t& seg, Point_t& origin, bool backwards) const
  {
    // Turn the trajectory into half-line that connect start -> end
    HalfLine_t lin1(trj.front(),trj.back()-trj.front());
    // Turn the segment into half-line
    HalfLine_t lin2(seg.Start(), seg.Dir());
    return _commonOrigin_(lin1, lin2, origin, backwards);
  }

  /// Common origin: Trajectory & Half Line. Keep track of origin
  double GeoAlgo::_commonOrigin_(const Trajectory_t& trj, const HalfLine_t& lin, Point_t& origin, bool backwards) const
  {
    // Turn the trajectory into half-line that connect start -> end
    HalfLine_t lin2(trj.front(),trj.back()-trj.front());
    return _commonOrigin_(lin, lin2, origin, backwards);
  }


  /// Bounding Sphere problem
  /// Real-Time Collision Analysis 4.3.5 (Pg. 100) - WelzlSphere
  Sphere_t GeoAlgo::_boundingSphere_(const std::vector<Point_t>& pts) const
  {

    // Remove any duplicate points
    std::vector<Point_t> copyPts = {pts[0]};
    for(size_t p1=0; p1 < pts.size(); p1++){
      // if an identical point does not already exist in copyPts -> then add
      bool found = false;
      for (size_t p2 =0; p2 < copyPts.size(); p2++)
	if (pts[p1] == copyPts[p2]) { found = true; break; }
      if (!found)
	copyPts.push_back(pts[p1]);
    }

    // if 4 or less points call appropriate constructor
    if (copyPts.size() < 5)
      return Sphere_t(copyPts);

    size_t npoints = copyPts.size();
    std::vector<Point_t> points4 = {copyPts[npoints-1],copyPts[npoints-2],copyPts[npoints-3],copyPts[npoints-4]};
    copyPts.pop_back();
    copyPts.pop_back();
    copyPts.pop_back();
    copyPts.pop_back();
    Sphere_t tmpSphere = Sphere(points4);
    return _RemainingPoints_(copyPts,tmpSphere);

    // too many points to call simple constructor! find minimally bounding sphere
    // compute sphere for first 4 points
    //Sphere_t tmpSphere(copyPts[0],copyPts[1],copyPts[2],copyPts[3]);
    //std::vector<Point_t> sosPoints;// = {pts[0]};
    //sosPoints.clear();
    //return _WelzlSphere_(copyPts,copyPts.size(),sosPoints);
  }

  Sphere_t GeoAlgo::_RemainingPoints_(std::vector<Point_t>& remaining,
				      const Sphere_t& thisSphere) const
  {


    // if no points lef -> done...return the current sphere
    if (remaining.size() == 0)
      return thisSphere;

    //std::cout << "Remaining points: " << remaining.size() << std::endl;

    auto const& lastPoint = remaining.back();

    // if this point is bounded by the already constructed sphere, continue
    if ( thisSphere.Contain(lastPoint) ){
      remaining.pop_back();
      return _RemainingPoints_(remaining,thisSphere);
    }
    // if not, need to adjust so that the new point is also bound
    // get distance from lastPoint and center
    double dist = lastPoint.Dist(thisSphere.Center());
    //std::cout << "point does not fit: " << lastPoint << std::endl;
    //std::cout << "distance: " << dist << std::endl;
    //std::cout << "center  : " << thisSphere.Center() << std::endl;
    //std::cout << "radius  : " << thisSphere.Radius() << std::endl;
    // the new center should be shifted in the direction
    // of the new point by half the difference between
    // the current radius and "dist"
    // direction in which to move:
    Vector_t dir = lastPoint-thisSphere.Center();
    dir.Normalize();
    // amount to move by
    double shift = (dist-thisSphere.Radius())/2.;
    if (shift < 0) { shift *= -1; }
    Point_t newCenter = thisSphere.Center() + dir*shift;
    double newRadius  = thisSphere.Radius()+shift;
    //std::cout << "new center: " << newCenter << std::endl;
    //std::cout << "new radius: " << newRadius << std::endl;
    Sphere_t newsphere(newCenter,newRadius);
    //if (newsphere.Contain(lastPoint)) { std::cout << "new point contained!" << std::endl; }
    //else { std::cout << "WRONG!" << std::endl; }
    remaining.pop_back();

    return _RemainingPoints_(remaining,newsphere);
  }

  Sphere_t GeoAlgo::_WelzlSphere_(const std::vector<Point_t>& pts,
				  int numPts,
				  std::vector<Point_t> sosPts) const
  {

    if (numPts == 0)
      return Sphere_t(sosPts);
    // choose last point in the input set as the one to test (if it fits in current sphere or not)
    int index = numPts-1;
    // recursively compute the smallest bounding sphere of the remaining points
    Sphere_t smallestSphere = _WelzlSphere_(pts, numPts-1, sosPts);
    // if the selected point lies inside this sphere, it is indeed the smallest
    if ( smallestSphere.Contain(pts[index]) )
      return smallestSphere;
    // otherwise, update the set-of-support to additionally contain the new point
    sosPts.push_back(pts[index]);
    return _WelzlSphere_(pts, numPts-1, sosPts);
  }

}
#endif
