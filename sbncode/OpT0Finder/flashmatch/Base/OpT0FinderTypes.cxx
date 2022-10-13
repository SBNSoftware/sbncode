#ifndef __OPT0FINDERTYPES_CXX__
#define __OPT0FINDERTYPES_CXX__

#include <iostream>
#include "OpT0FinderTypes.h"
#include <math.h> 

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

  double QCluster_t::pitch() const
  {
    double pitch=0.3;
    // find first track point
    auto init_pt = (*this)[0];
    for (auto const& pt : (*this)){
      if (pt.trk!=1) continue;
      else{ 
        init_pt = pt;
        break;
      }
    }
    double min_x = init_pt.x; 
    double max_x = init_pt.x;
    double min_z = init_pt.z;
    double max_z = init_pt.z;
    for (auto const& pt : (*this)){
      // std::cout << pt.trk << std::endl;
      if(pt.trk!= 1) continue;// if the point does not belong to a track, skip
      if(pt.x > max_x) max_x = pt.x;
      if(pt.x < min_x) min_x = pt.x;
      if(pt.z > max_z) max_z = pt.z;
      if(pt.z < min_z) min_z = pt.z;
    }
    // std::cout << "x: (min, max) // z: (min, max): (" << min_x << ", " << max_x << "), (" << min_z << ", " << max_z << ")" << std::endl;
    double gamma = atan(abs(max_x - min_x) / abs(max_z - min_z));
    // std::cout << "gamma: " << gamma << std::endl;
    pitch = 0.3/cos(gamma); 
    // std::cout << "pitch: " << pitch << std::endl;
    return pitch;
  }

  /// streamer override
  std::ostream& operator << (std::ostream& out, const flashmatch::QCluster_t& obj) {
    out << "QCluster_t " << obj.size() << " points length=" << obj.length() << " qsum=" << obj.sum() << std::endl;
    return out;
  }




}

#endif
