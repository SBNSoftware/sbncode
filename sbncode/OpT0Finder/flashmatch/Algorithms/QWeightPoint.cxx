#ifndef OPT0FINDER_QWEIGHTPOINT_CXX
#define OPT0FINDER_QWEIGHTPOINT_CXX

#include "QWeightPoint.h"
#include <cmath>
#include <sstream>
#include <numeric>
namespace flashmatch {

  static QWeightPointFactory __global_QWeightPointFactory__;

  QWeightPoint::QWeightPoint(const std::string name)
    : BaseFlashMatch(name)
    , _x_step_size ( 0.5   )
    , _zdiff_max   ( 50*50 )
  {}

  void QWeightPoint::_Configure_(const Config_t &pset)
  {
    _x_step_size = pset.get<double>("XStepSize");
    _zdiff_max   = pset.get<double>("ZDiffMax" );
    _zdiff_max *= _zdiff_max;
  }

  void QWeightPoint::Match(const QCluster_t& pt_v, const Flash_t& flash, FlashMatch_t& match)
  {

    if(_vis_array.pe_v.empty())
      _vis_array.pe_v.resize(DetectorSpecs::GetME().NOpDets());

    // Prepare the return values (Mostly QWeightPoint)
    if(pt_v.empty()){
      std::cout<<"Not enough points!"<<std::endl;
      return;
    }

    _tpc_qcluster.resize(pt_v.size());

    // Get min & max x value
    double x_max = 0;
    double x_min = 1e12;

    for(auto const& pt : pt_v) {
      if(pt.x > x_max) x_max = pt.x;
      if(pt.x < x_min) x_min = pt.x;
    }

    double min_dz = 1e9;
    for(double x_offset=0; x_offset<(256.35-(x_max-x_min)); x_offset+=_x_step_size) {

      // Create QCluster_t with this offset

      for(size_t i=0; i<_tpc_qcluster.size(); ++i) {
	_tpc_qcluster[i].x = pt_v[i].x + x_offset - x_min;
	_tpc_qcluster[i].y = pt_v[i].y;
	_tpc_qcluster[i].z = pt_v[i].z;
	_tpc_qcluster[i].q = pt_v[i].q;
      }

      FillEstimate(_tpc_qcluster,_vis_array);

      // Calculate amplitudes corresponding to max opdet amplitudes
      double vis_pe_sum = _vis_array.TotalPE();

      double weighted_z = 0;
      for(size_t pmt_index=0; pmt_index<DetectorSpecs::GetME().NOpDets(); ++pmt_index) {

	if(_vis_array.pe_v[pmt_index]<0) continue;
	weighted_z += DetectorSpecs::GetME().PMTPosition(pmt_index)[2] * _vis_array.pe_v[pmt_index] / vis_pe_sum;

      }

      double dz = std::fabs(weighted_z - flash.z);

      if(dz < min_dz) {

	min_dz = dz;

	match.score = 1./min_dz;
	match.tpc_point.x = match.tpc_point.y = 0;
	match.tpc_point.q = vis_pe_sum;

	match.tpc_point.x = x_offset;

	for(size_t pmt_index=0; pmt_index<DetectorSpecs::GetME().NOpDets(); ++pmt_index) {
	  if(_vis_array.pe_v[pmt_index]<0) continue;
	  match.tpc_point.y += DetectorSpecs::GetME().PMTPosition(pmt_index)[1] * _vis_array.pe_v[pmt_index] / vis_pe_sum;
	}

	match.tpc_point.z = weighted_z;
      }
    }

    match.hypothesis.clear();

    FLASH_INFO() << "Best match Hypothesis: "
		 << match.tpc_point.x << " : "
		 << match.tpc_point.y << " : "
		 << match.tpc_point.z << " ... min dist : " << min_dz
		 << std::endl;

    // If min-diff is bigger than assigned max, return default match (score<0)
    if( min_dz > _zdiff_max ) {
      match.tpc_point.x = match.tpc_point.y = match.tpc_point.z = -1;
      match.tpc_point.q = -1;
      match.score = -1;
      return;
    }

    match.hypothesis = _vis_array.pe_v;

  }


}
#endif
