#ifndef LightPath_CXX
#define LightPath_CXX

#include "LightPath.h"

namespace flashmatch {

  static LightPathFactory __global_LightPathFactory__;

  LightPath::LightPath(const std::string name)
    : BaseAlgorithm(kCustomAlgo, name)
    , _gap         ( 0.5    )
    , _light_yield ( 40000. )
    , _dEdxMIP     ( 2.07   ) //1.42[Mev*cm^2*g]*1.4[g/cm^3]=2.004MeV/cm
  {}

  void LightPath::_Configure_(const Config_t &pset)
  {
    _gap          = pset.get< double > ( "SegmentSize" );
    _light_yield  = pset.get< double > ( "LightYield"  );
    _dEdxMIP      = pset.get< double > ( "MIPdEdx"     );
  }

  void LightPath::MakeQCluster(const ::geoalgo::Vector& pt_1,
			       const ::geoalgo::Vector& pt_2,
			       QCluster_t& Q_cluster,
			       double dedx) const {

    if(dedx < 0) dedx = _dEdxMIP;

    double dist = pt_1.Dist(pt_2);
    QPoint_t q_pt;
    FLASH_INFO() << "Filling points between (" << pt_1[0] << "," << pt_1[1] << "," << pt_1[2] << ")"
		 << " => (" << pt_2[0] << "," << pt_2[1] << "," << pt_2[2] << ") ... dist="<<dist<<std::endl;
    if (dist <= _gap) {
      ::geoalgo::Vector mid_pt((pt_1 + pt_2) / 2.);
      q_pt.x = mid_pt[0];
      q_pt.y = mid_pt[1];
      q_pt.z = mid_pt[2];
      q_pt.q = _light_yield * _dEdxMIP * dist;
      q_pt.trk = -1;
      FLASH_DEBUG() << "Smaller than gap threshold (" << _gap << ")" << std::endl
		   << "Traj pt (" << q_pt.x << "," << q_pt.y << "," << q_pt.z << ") q=" << q_pt.q << std::endl;
      Q_cluster.emplace_back(q_pt);
      return;
    }

    int num_div = int(dist / _gap);

    ::geoalgo::Vector direct = (pt_1 - pt_2).Dir();

    Q_cluster.reserve(Q_cluster.size() + num_div);

    for (int div_index = 0; div_index < num_div + 1; div_index++) {
      if (div_index < num_div) {
        auto const mid_pt = pt_2 + direct * (_gap * div_index + _gap / 2.);
        q_pt.x = mid_pt[0] ;
        q_pt.y = mid_pt[1];
        q_pt.z = mid_pt[2];
        q_pt.q = _light_yield * _dEdxMIP * _gap;
	FLASH_DEBUG() << "Traj pt (" << q_pt.x << "," << q_pt.y << "," << q_pt.z << ") q=" << q_pt.q << std::endl;
        Q_cluster.emplace_back(q_pt);
      }
      else {
        double weight = (dist - int(dist / _gap) * _gap);
        auto const mid_pt = pt_2 + direct * (_gap * div_index + weight / 2.);
        q_pt.x = mid_pt[0] ;
        q_pt.y = mid_pt[1];
        q_pt.z = mid_pt[2];
        q_pt.q = _light_yield * _dEdxMIP * weight;
	FLASH_DEBUG() << "Traj pt (" << q_pt.x << "," << q_pt.y << "," << q_pt.z << ") q=" << q_pt.q << std::endl;
        Q_cluster.emplace_back(q_pt);
      }//Last segment less than gap
    }
  }

  QCluster_t LightPath::MakeQCluster(const ::geoalgo::Trajectory& trj) const {

    QCluster_t result;
    result.clear();

    for (size_t i = 0; i < trj.size() - 1; i++) {
      auto const& this_loc(trj[i]);
      auto const& last_loc(trj[i + 1]);
      LightPath::MakeQCluster(this_loc, last_loc, result);
    }

    // Trimming Q_cluster
    auto const& bbox = DetectorSpecs::GetME().ActiveVolume();
    double _vol_xmax = bbox.Max()[0];
    double _vol_ymax = bbox.Max()[1];
    double _vol_zmax = bbox.Max()[2];

    double _vol_xmin = bbox.Min()[0];
    double _vol_ymin = bbox.Min()[1];
    double _vol_zmin = bbox.Min()[2];
    FLASH_INFO() << result << std::endl;
    QCluster_t final_result;
    final_result.clear();
    for (size_t idx = 0; idx < result.size(); ++idx) {
        auto pt = result[idx];
        if (pt.x >= _vol_xmin && pt.x <= _vol_xmax && pt.y >= _vol_ymin && pt.y <= _vol_ymax && pt.z >= _vol_zmin && pt.z <= _vol_zmax) {
            final_result.push_back(pt);
        }
    }
    FLASH_INFO() << final_result << std::endl;

    return final_result;
  }

}


#endif
