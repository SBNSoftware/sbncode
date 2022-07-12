#ifndef OPT0FINDER_AnodeMatch_CXX
#define OPT0FINDER_AnodeMatch_CXX

#include "AnodeMatch.h"
#include <cmath>
#include <sstream>
#include <numeric>
namespace flashmatch {

  static AnodeMatchFactory __global_AnodeMatchFactory__;

  AnodeMatch::AnodeMatch(const std::string name)
    : BaseTouchMatch(name)
  {}

  void AnodeMatch::_Configure_(const Config_t &pset)
  {
    _time_window  = pset.get<double>("TimeWindow");
    _pe_threshold = pset.get<double>("PEThreshold",-1);
    _space_range  = pset.get<double>("DistanceToPMTs",50);
  }

  void AnodeMatch::Match(const QCluster_t& pt_v, const Flash_t& flash, FlashMatch_t& match)
  {

    FLASH_INFO() << "Starting" << std::endl;

    double detector_xmin = DetectorSpecs::GetME().ActiveVolume().Min()[0];
    double detector_xmax = DetectorSpecs::GetME().ActiveVolume().Max()[0];

    QPoint_t qpt_min, qpt_max;
    pt_v.min_x(qpt_min);
    pt_v.max_x(qpt_max);

    FLASH_INFO() << "TPC number of points " << pt_v.size() << std::endl;
    FLASH_INFO() << "TPC min pt: " << qpt_min.x << " " << qpt_min.y << " " << qpt_min.z << std::endl;
    FLASH_INFO() << "TPC max pt: " << qpt_max.x << " " << qpt_max.y << " " << qpt_max.z << std::endl;

    // Calculate the hypothesized flash time assuming the track is touching the TPC min-x boundary
    double time_tpc0 = (qpt_min.x - detector_xmin) / DetectorSpecs::GetME().DriftVelocity();

    // Calculate the hypothesized flash time assuming the track is touching the TPC max-x boundary
    double time_tpc1 = (detector_xmax - qpt_max.x) / DetectorSpecs::GetME().DriftVelocity();

    FLASH_INFO() << "Timing if tuching TPC0 : " << time_tpc0  << std::endl;
    FLASH_INFO() << "Timing if tuching TPC1 : " << time_tpc1  << std::endl;
    FLASH_INFO() << "Flash timing           : " << flash.time << std::endl;

    // Check if each time can be consistent with this flash
    double time_proximity_tpc0 = std::fabs(time_tpc0 - flash.time);
    double time_proximity_tpc1 = std::fabs(flash.time - time_tpc1);

    FLASH_INFO() << "TPC0 time proximity : " << time_tpc0 - flash.time << std::endl;
    FLASH_INFO() << "TPC1 time proximity : " << time_tpc1 - flash.time << std::endl;
    FLASH_INFO() << "Time window         : " << _time_window << std::endl;
    // Check TPC0 boundary
    bool touching_tpc0 = false;
    double touching_dist0 = kINVALID_DOUBLE;
    if(time_proximity_tpc0 < _time_window) {
      FLASH_INFO() << "Time compatible with TPC 0..." << std::endl;
      // time is compatible. next, check the PE amount in the nearby PMT
      ::geoalgo::Point_t tpc_pt(detector_xmin, qpt_min.y, qpt_min.z);
      double hist_min_dist = kINVALID_DOUBLE;
      double hist_pe = -1;
      int hist_closest_pmt = -1;
      for(size_t i=0; i<DetectorSpecs::GetME().NOpDets(); ++i) {
        auto const& pmt_pos = DetectorSpecs::GetME().PMTPosition(i);
        double dist = tpc_pt.Dist(pmt_pos);
        if(dist < hist_min_dist) {
          hist_min_dist = dist;
          hist_pe = flash.pe_v[i];
          hist_closest_pmt = i;
        }
        if(dist > _space_range) continue;
        if(flash.pe_v[i] < _pe_threshold) continue;
        // PMT hit is found that is consistent with a touch match
        touching_tpc0 = true;
        touching_dist0 = std::min(dist,touching_dist0);
        FLASH_INFO() << "TOUCH FOUND ... TPC 0 PMT " << i << " dist " << sqrt(dist) << " ... " << flash.pe_v[i] << " PE" <<std::endl;
      }
      FLASH_INFO() << "... in TPC 0, closest PMT " << hist_closest_pmt << " dist. " << hist_min_dist << " PE " << hist_pe << std::endl;
    }

    // Check TPC1 boundary
    bool touching_tpc1 = false;
    double touching_dist1 = kINVALID_DOUBLE;
    if(time_proximity_tpc1 < _time_window) {
      FLASH_INFO() << "Time compatible with TPC 1..." << std::endl;
      // time is compatible. next, check the PE amount in the nearby PMT
      ::geoalgo::Point_t tpc_pt(detector_xmax, qpt_max.y, qpt_max.z);
      double hist_min_dist = kINVALID_DOUBLE;
      double hist_pe = -1;
      int hist_closest_pmt = -1;
      for(size_t i=0; i<DetectorSpecs::GetME().NOpDets(); ++i) {
        auto const& pmt_pos = DetectorSpecs::GetME().PMTPosition(i);
        double dist = tpc_pt.Dist(pmt_pos);
        if(dist < hist_min_dist) {
          hist_min_dist = dist;
          hist_pe = flash.pe_v[i];
          hist_closest_pmt = i;
        }
        if(dist > _space_range) continue;
        if(flash.pe_v[i] < _pe_threshold) continue;
        // PMT hit is found that is consistent with a touch match
        touching_tpc1 = true;
        touching_dist1 = std::min(dist,touching_dist1);
        FLASH_INFO() << "TOUCH FOUND ... TPC 1 PMT " << i << " dist " << sqrt(dist) << " ... " << flash.pe_v[i] << " PE" <<std::endl;
      }
      FLASH_INFO() << "... in TPC 1, closest PMT " << hist_closest_pmt << " dist. " << hist_min_dist << " PE " << hist_pe << std::endl;
    }

    if(touching_tpc0 || touching_tpc1) {
      match.touch_match = flashmatch::kAnodeCrossing;
      match.touch_score = std::min(time_proximity_tpc0,time_proximity_tpc1);
      match.touch_point = (touching_dist0 < touching_dist1 ? qpt_min : qpt_max);
    }


  }


}
#endif
