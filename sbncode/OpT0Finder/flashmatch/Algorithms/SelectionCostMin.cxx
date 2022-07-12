#ifndef OPT0FINDER_SelectionCostMin_CXX
#define OPT0FINDER_SelectionCostMin_CXX

#include "SelectionCostMin.h"

namespace flashmatch {
  
  static SelectionCostMinFactory __global_SelectionCostMinFactory__;

  SelectionCostMin::SelectionCostMin(const std::string name)
    : BaseMatchSelection(name)
  {}

  void SelectionCostMin::_Configure_(const Config_t &pset)
  {
  }  


  std::vector<FlashMatch_t> 
  SelectionCostMin::Select(const std::vector<std::vector<FlashMatch_t> >& match_data)
  {
    std::vector<FlashMatch_t> result;
    return result;
  }

}

#endif
