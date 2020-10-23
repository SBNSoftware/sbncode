#ifndef OPT0FINDER_COMMONAMPS_CXX
#define OPT0FINDER_COMMONAMPS_CXX

#include "CommonAmps.h"
#include "flashmatch/Base/OpT0FinderException.h"
#include <cmath>
#include <sstream>
#include <numeric>
namespace flashmatch {

  static CommonAmpsFactory __global_CommonAmpsFactory__;

  CommonAmps::CommonAmps(const std::string name)
    : BaseFlashMatch(name)
  {
    _percent = 0.5;
    _score   = 0.8;    
  }

  void CommonAmps::_Configure_(const Config_t &pset)
  {
    _percent = pset.get<double>("QFracThreshold");
    _score   = pset.get<double>("ScoreThreshold");
    _x_step_size = pset.get<double>("XStepSize");
  }
  
  FlashMatch_t CommonAmps::Match(const QCluster_t& pt_v, const Flash_t& flash)
  {
    
    double integral_op  = std::accumulate(std::begin(flash.pe_v),
					  std::end(flash.pe_v),
					  0.0);
    double maxRatio 	= -1;
    //double maxX         = 0;
    
    if(_vis_array.pe_v.empty())
      _vis_array.pe_v.resize(OpDetXArray().size());

    // Create multimap to hold the largest amplitudes in the first slots of map
    // Normalize each PE bin to the total # of PEs-- 1/x to put highest amps in front
    std::multimap<double,int> ampToOpDet ;

    for(size_t k=0; k<32; k++)
      ampToOpDet.emplace(1./(flash.pe_v[k]/integral_op),k);
    

    // Now that opdet hits are normalized and ordered, store indices for PMTs which 
    // had the greatest hit amplitudes up to _percent of total PMT hits
    double opAmpTotal = 0;
    std::vector<int> ids ;
    ids.resize(0);

    for( auto const & e : ampToOpDet ){
      opAmpTotal += 1/e.first ;
      ids.push_back(e.second) ;
      if( opAmpTotal > _percent )
	break;
    }

    // Prepare the return values (Mostly QWeightPoint)
    FlashMatch_t f;
    if(pt_v.empty()){
      std::cout<<"Not enough points!"<<std::endl;
      return f;
    }
    
    _tpc_qcluster.resize(pt_v.size());

    // Get min & max x value
    double x_max = 0;
    double x_min = 1e12;

    for(auto const& pt : pt_v) {
      if(pt.x > x_max) x_max = pt.x;
      if(pt.x < x_min) x_min = pt.x;
    }

    for(double x_offset=0;
	x_offset<(250.-(x_max-x_min));
	x_offset+=_x_step_size) {
      
      // Create QCluster_t with this offset

      for(size_t i=0; i<_tpc_qcluster.size(); ++i) {
	_tpc_qcluster[i].x = pt_v[i].x + x_offset - x_min;
	_tpc_qcluster[i].y = pt_v[i].y;
	_tpc_qcluster[i].z = pt_v[i].z;
	_tpc_qcluster[i].q = pt_v[i].q;
      }
            
      FillEstimate(_tpc_qcluster,_vis_array);

      // Calculate amplitudes corresponding to max opdet amplitudes
      double visAmpTotal = 0;
      double vis_pe_sum = _vis_array.TotalPE();

      for(size_t i=0; i<ids.size(); i++) {
	if(_vis_array.pe_v[ids[i]]<0) continue;
	visAmpTotal += _vis_array.pe_v[ids[i]] / vis_pe_sum ;
      }
      
      double ratio = 0;
      if(opAmpTotal > visAmpTotal )
	ratio = visAmpTotal/opAmpTotal ;
      else
	ratio = opAmpTotal/visAmpTotal ;
      
      if(ratio > maxRatio) {
	maxRatio = ratio;
	//maxX     = x_offset;

	f.score = ratio;
	f.tpc_point.x = f.tpc_point.y = f.tpc_point.z = 0;
	f.tpc_point.q = vis_pe_sum;

	for(size_t pmt_index=0; pmt_index<NOpDets(); ++pmt_index) {
	  double pe = _vis_array.pe_v[pmt_index];
	  if(pe<0) continue;
	  f.tpc_point.x += OpDetX(pmt_index) * pe / vis_pe_sum;
	  f.tpc_point.y += OpDetY(pmt_index) * pe / vis_pe_sum;
	  f.tpc_point.z += OpDetZ(pmt_index) * pe / vis_pe_sum;
	}
      }
    }
  
    // If min-diff is bigger than assigned max, return default match (score<0)
    if( maxRatio < _score ) {

      f.tpc_point.x = f.tpc_point.y = f.tpc_point.z = -1;
      f.tpc_point.q = -1;
      f.score = -1;
      return f;
    }

    f.hypothesis = _vis_array.pe_v;
    return f;
  }
}
#endif
