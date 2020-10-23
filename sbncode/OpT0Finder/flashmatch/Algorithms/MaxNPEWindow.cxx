#ifndef OPT0FINDER_MAXNPEWINDOW_CXX
#define OPT0FINDER_MAXNPEWINDOW_CXX

#include "MaxNPEWindow.h"
#include "TimeRange.h"
#include <map>
#include <numeric>
//#include <functional>
namespace flashmatch {

  static MaxNPEWindowFactory __global_MaxNPEWindowFactory__;

  MaxNPEWindow::MaxNPEWindow(const std::string name)
    : BaseFlashFilter(name)
    , _time_lower_bound(-0.1) // Default T1 in micro-second
    , _time_upper_bound(8.)   // Default T2 in micro-second
    , _npe_threshold(10)       // Default # p.e. as a threshold
  {}

  void MaxNPEWindow::_Configure_(const Config_t &pset)
  {
    _time_lower_bound = pset.get<double>("TimeLowerBound");
    _time_upper_bound = pset.get<double>("TimeUpperBound");
    _npe_threshold    = pset.get<double>("NPEThreshold"  );
  }

  IDArray_t MaxNPEWindow::Filter(const FlashArray_t& flash_v)
  {
    // order flash via npe
    std::multimap<double,ID_t> npe_map;
    // Loop over flash array
    for(size_t index=0; index<flash_v.size(); ++index) {

      auto const& flash = flash_v[index]; // Retrieve this flash

      double npe = std::accumulate(flash.pe_v.begin(),flash.pe_v.end(),0.0); // Sum p.e.

      if(npe < _npe_threshold) continue; // Ignore if below threshold
      
      npe_map.emplace( 1./npe, index ); // Use inverse of npe to order high=>low

    }

    // Prepare a return flashmatch::IDArray_t object
    IDArray_t ordered_index_v;
    ordered_index_v.reserve(npe_map.size()); // reserve for fast push_back

    TimeRangeSet trs; // Use TimeRange utility to detect overlap windows
    for(auto const& npe_index : npe_map) {

      auto const& flash = flash_v[npe_index.second]; // retrieve flash

      if(trs.Overlap(flash.time)) continue; // if within defined set of windows, ignore

      ordered_index_v.push_back(npe_index.second); // candidate: add to a list

      // update TimeRangeSet to keep track of windows.
      trs.Insert( TimeRange(flash.time + _time_lower_bound,
			    flash.time + _time_upper_bound) );
    }

    return ordered_index_v;
  }


}

#endif
