#ifndef OPT0FINDER_TIMECOMPATMATCH_CXX
#define OPT0FINDER_TIMECOMPATMATCH_CXX

#include "TimeCompatMatch.h"
#include <cmath>
#include <sstream>

namespace flashmatch {

  static TimeCompatMatchFactory __global_TimeCompatMatchFactory__;

  TimeCompatMatch::TimeCompatMatch(const std::string name)
    : BaseProhibitAlgo(name)
  {}

  void TimeCompatMatch::_Configure_(const Config_t &pset)
  {
    _time_shift  = pset.get<double>("BeamTimeShift",0);
    _time_window = pset.get<double>("TouchingTrackWindow",5);
  }

  bool TimeCompatMatch::MatchCompatible(const QCluster_t& clus, const Flash_t& flash)
  {

    if(clus.empty()) return false; 

    // get time of flash
    auto flash_time = flash.time - _time_shift;

    // get time of cluster by looking at the range of x-positions
    double clus_x_min = kINVALID_DOUBLE;
    double clus_x_max = -1 * kINVALID_DOUBLE;
    for (auto const& pt : clus){
      if (pt.x > clus_x_max) { clus_x_max = pt.x; }
      if (pt.x < clus_x_min) { clus_x_min = pt.x; }
    }

    // Detector boundary info
    double xmax = DetectorSpecs::GetME().ActiveVolume().Max()[0];
    double xmin = DetectorSpecs::GetME().ActiveVolume().Min()[0];

    // Assume this is the right flash, check if a trajectory is within the active volume
    // One assumes tpc0 and the other tpc1 
    double reco_x_tpc0 = clus_x_min - flash_time * DetectorSpecs::GetME().DriftVelocity();
    double reco_x_tpc1 = clus_x_max + flash_time * DetectorSpecs::GetME().DriftVelocity();

    double distance_window = _time_window * DetectorSpecs::GetME().DriftVelocity();

    bool incompatible_tpc0 = (xmin - reco_x_tpc0) > distance_window;
    bool incompatible_tpc1 = (reco_x_tpc1 - xmax) > distance_window;

    FLASH_INFO() << "Inspecting..." << std::endl
      << "Detector X span : " << xmin << " => " << xmax << std::endl
      << "TPC pts X span  : " << clus_x_min << " => " << clus_x_max << " ... " << clus.size() << " points" << std::endl
      << "Flash time      : " << flash_time << " (shifted by " << _time_shift << ")" << std::endl
      << "Hypothesis X pos: " << reco_x_tpc0 << " => " << reco_x_tpc1 << std::endl
      << "From TPC-0 edge : " << (xmin - reco_x_tpc0) << " ... incompatible? " << incompatible_tpc0 << std::endl
      << "From TPC-1 edge : " << (reco_x_tpc1 - xmax) << " ... incompatible? " << incompatible_tpc1 << std::endl;

    if ( incompatible_tpc0 && incompatible_tpc1 ) return false;
    return true;

  }


}
#endif
