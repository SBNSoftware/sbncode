#ifndef OPT0FINDER_NPTFILTER_CXX
#define OPT0FINDER_NPTFILTER_CXX

#include "NPtFilter.h"

namespace flashmatch {
  
  static NPtFilterFactory __global_NPtFilterFactory__;

  NPtFilter::NPtFilter(const std::string name)
    : BaseTPCFilter(name)
  { _min_num_pt = 2; }

  void NPtFilter::_Configure_(const Config_t &pset)
  {
    _min_num_pt = pset.get<double>("MinNumPoint");
  }
  
  IDArray_t NPtFilter::Filter(const QClusterArray_t& tpc_obj_v) {

    IDArray_t res;

    // Loop over QCluster_t array
    for(ID_t id=0; id<tpc_obj_v.size(); ++id) {

      auto const& tpc_obj = tpc_obj_v[id]; // retrieve

      // if more # of QPoint_t than threshold, accept
      if(tpc_obj.size() >= _min_num_pt) res.push_back(id); 

    }

    return res;
  }
  

}

#endif
