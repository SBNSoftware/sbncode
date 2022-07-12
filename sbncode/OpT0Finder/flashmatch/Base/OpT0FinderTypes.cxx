#ifndef __OPT0FINDERTYPES_CXX__
#define __OPT0FINDERTYPES_CXX__

#include <iostream>
#include "OpT0FinderTypes.h"

namespace flashmatch {

  void QCluster_t::max_x(flashmatch::QPoint_t& qpt) const 
  {
    double max_value = -kINVALID_DOUBLE;
    for(auto const& p : (*this)) {
      if(p.x < max_value) continue;
      max_value = p.x;
      qpt = p;      
    }
  }

  void QCluster_t::max_x(geoalgo::Point_t& pt) const 
  {
    QPoint_t qpt;
    this->max_x(qpt);
    pt[0] = qpt.x;
    pt[1] = qpt.y;
    pt[2] = qpt.z;
  }

  void QCluster_t::min_x(flashmatch::QPoint_t& qpt) const 
  {
    double min_value = kINVALID_DOUBLE;
    for(auto const& p : (*this)) {
      if(p.x > min_value) continue;
      min_value = p.x;
      qpt = p;      
    }
  }

  void QCluster_t::min_x(geoalgo::Point_t& pt) const 
  {
    QPoint_t qpt;
    this->min_x(qpt);
    pt[0] = qpt.x;
    pt[1] = qpt.y;
    pt[2] = qpt.z;
  }

  double QCluster_t::sum() const
  { double sum=0; for(auto const& pt : (*this)) sum += pt.q; return sum; }

  double QCluster_t::length() const
  {
    double len=0.;
    for(size_t idx=1; idx<this->size(); ++idx) {
      auto const& pt0 = (*this)[idx-1];
      auto const& pt1 = (*this)[idx];
      len += sqrt(pow(pt0.x - pt1.x,2)+pow(pt0.y - pt1.y,2)+pow(pt0.z - pt1.z,2));
    }
    return len;
  }

  void QCluster_t::drop(double x_min, double x_max)
  {
    QCluster_t another;
    another.idx = this->idx;
    another.time = this->time;
    another.time_true = this->time_true;
    another.min_x_true = this->min_x_true;
    another.reserve(this->size());
    for(auto const& pt : (*this)) {
      if(pt.x < x_min) continue;
      if(pt.x > x_max) continue;
      another.push_back(pt);
    }
    (*this) = another;
  }

  void QCluster_t::drop(double x_min, double y_min, double z_min,
                        double x_max, double y_max, double z_max) {
    QCluster_t another;
    another.idx = this->idx;
    another.time = this->time;
    another.time_true = this->time_true;
    another.min_x_true = this->min_x_true;
    another.reserve(this->size());
    for(auto const& pt : (*this)) {
      if(pt.x < x_min) continue;
      if(pt.x > x_max) continue;
      if(pt.y < y_min) continue;
      if(pt.y > y_max) continue;
      if(pt.z < z_min) continue;
      if(pt.z > z_max) continue;
      another.push_back(pt);
    }
    (*this) = another;
  }

  /// streamer override
  std::ostream& operator << (std::ostream& out, const flashmatch::QCluster_t& obj) {
    out << "QCluster_t " << obj.size() << " points length=" << obj.length() << " qsum=" << obj.sum() << std::endl;
    return out;
  }




}

#endif
