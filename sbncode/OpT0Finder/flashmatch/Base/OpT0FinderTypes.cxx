#ifndef __OPT0FINDERTYPES_CXX__
#define __OPT0FINDERTYPES_CXX__

#include <iostream>
#include "OpT0FinderTypes.h"

namespace flashmatch {
  
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
    another.reserve(this->size());
    for(auto const& pt : (*this)) {
      if(pt.x < x_min) continue;
      if(pt.x > x_max) continue;
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
