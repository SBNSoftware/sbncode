#ifndef OPT0FINDER_FILTERARRAY_CXX
#define OPT0FINDER_FILTERARRAY_CXX

#include "FilterArray.h"

namespace flashana {

  FilterArray::FilterArray() : BaseTPCFilter()
  { }
  
  IDArray_t FilterArray::Filter(const QClusterArray_t& tpc_obj_v) {

    // list of indices to be kept
    IDArray_t res;

    // temporary tpc_obj vector
    QClusterArray_t tmp_tpc_obj_v;
    for(ID_t id=0; id<tpc_obj_v.size(); ++id) {
      tmp_tpc_obj_v.push_back( tpc_obj_v[id] );
      res.push_back(id);
    }

    // apply each algorithm sequentially
    for (size_t i=0; i < _filter_v.size(); i++){

      auto indices = _filter_v[i]->Filter(tmp_tpc_obj_v);

      // make a copy of the result indices and the tpc_obj vector
      auto res_cpy = res;
      auto tpc_obj_cpy_v = tmp_tpc_obj_v;
      res.clear();
      tmp_tpc_obj_v.clear();

      for (auto& idx : indices){
	res.push_back(res_cpy[idx]);
	tmp_tpc_obj_v.push_back(tpc_obj_cpy_v[idx]);
      }

    }// for all filter algorithms to apply

    return res;
  }
  

}

#endif
