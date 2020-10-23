#ifndef OPT0FINDER_TIMECOMPATMATCH_CXX
#define OPT0FINDER_TIMECOMPATMATCH_CXX

#include "TimeCompatMatch.h"

#ifndef USING_LARSOFT
#define USING_LARSOFT 1
#endif

namespace flashmatch {

  static TimeCompatMatchFactory __global_TimeCompatMatchFactory__;

  TimeCompatMatch::TimeCompatMatch(const std::string name)
    : BaseProhibitAlgo(name)
  {}

  void TimeCompatMatch::_Configure_(const Config_t &pset)
  {
    _time_buffer = pset.get<double>("TimeBuffer");
  }

  bool TimeCompatMatch::MatchCompatible(const QCluster_t& clus, const Flash_t& flash)
  {
    if(clus.empty()) {
      FLASH_INFO() << "QCluster_t is empty." << std::endl;
      return false;
    }

    // get time of flash
    auto flash_time = flash.time;

    // get time of cluster by looking at the range of x-positions
    // FIXME: 1036?
    double clus_x_min =  1036.; // cm
    double clus_x_max = -1036.;    // cm
    for (auto const& pt : clus){
      if (pt.x > clus_x_max) { clus_x_max = pt.x; }
      if (pt.x < clus_x_min) { clus_x_min = pt.x; }
    }

    FLASH_INFO() << "Cluster x min: " << clus_x_min << ", x max: " << clus_x_max << std::endl;

    // Earliest flash time => assume clus_x_max is @ detector X-max boundary
    #if USING_LARSOFT == 1
    double xmax = DetectorSpecs::GetME().ActiveVolume().Max()[0];
    double clus_t_min = (clus_x_max - xmax) / DetectorSpecs::GetME().DriftVelocity();
    double clus_t_max = clus_x_min / DetectorSpecs::GetME().DriftVelocity();
    #else
    double xmax = DetectorSpecs::GetME().ActiveVolume().Max()[0];
    double clus_t_min = (clus_x_max - xmax) / DetectorSpecs::GetME().DriftVelocity();
    double clus_t_max = clus_x_min / DetectorSpecs::GetME().DriftVelocity();
    #endif
    FLASH_INFO() << "Cluster xmax: " << xmax << ", clus_t_min: " << clus_t_min
                                             << ", clus_t_max: " << clus_t_max
                                             << ", flash_time: " << flash_time << std::endl;

    /*
    std::cout<< "Inspecting TPC object @ " << clus.time << std::endl;
    std::cout<< "xmin = " << clus_x_min << " ... xmax = " << clus_x_max << std::endl;
    std::cout<< "tmin = " << clus_t_min << " ... tmax = " << clus_t_max << std::endl;
    std::cout<< "Flash time @ " << flash_time << std::endl;
    */
    return ((clus_t_min - _time_buffer) < flash_time && flash_time < (clus_t_max + _time_buffer));

  }


}
#endif
