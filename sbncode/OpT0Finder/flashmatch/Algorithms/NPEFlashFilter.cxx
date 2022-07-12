#ifndef OPT0FINDER_NPEFLASHFILTER_CXX
#define OPT0FINDER_NPEFLASHFILTER_CXX

#include "NPEFlashFilter.h"
#include "TimeRange.h"
#include <map>
#include <numeric>
//#include <functional>
namespace flashmatch {

  static NPEFlashFilterFactory __global_NPEFlashFilterFactory__;

  NPEFlashFilter::NPEFlashFilter(const std::string name)
    : BaseFlashFilter(name)
    , _npe_threshold(10)       // Default # p.e. as a threshold
  {}

  void NPEFlashFilter::_Configure_(const Config_t &pset)
  {
    _npe_threshold    = pset.get<double>("NPEThreshold"  );
  }

  IDArray_t NPEFlashFilter::Filter(const FlashArray_t& flash_v)
  {
    // Prepare a return flashmatch::IDArray_t object
    IDArray_t result;
    
    // Loop over flash array
    for(size_t index=0; index<flash_v.size(); ++index) {

      auto const& flash = flash_v[index]; // Retrieve this flash

      double npe = std::accumulate(flash.pe_v.begin(),flash.pe_v.end(),0.0); // Sum p.e.

      if(npe < _npe_threshold) continue; // Ignore if below threshold

      result.push_back(index);

    }

    return result;
  }


}

#endif
