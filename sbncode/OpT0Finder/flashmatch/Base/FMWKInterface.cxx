#ifndef __OPT0FINDERFMWKINTERFACE_CXX__
#define __OPT0FINDERFMWKINTERFACE_CXX__

#include "FMWKInterface.h"
#include <assert.h>

namespace flashmatch{
  DetectorSpecs* DetectorSpecs::_me = nullptr;
}

namespace flashmatch{

  void DetectorSpecs::DumpInfo() const {
    FLASH_NORMAL() << "Information dump..." << std::endl << std::endl 
    << "  Active    volume : " << _bbox << std::endl << std::endl 
    << "  PhotonLib volume : " << _photon_bbox << std::endl << std::endl 
    << "  Drift velocity   : " << DriftVelocity() << std::endl 
    << "  LightYield       : " << LightYield() << std::endl 
    << "  MIP dE/dX        : " << _MIPdEdx << std::endl 
    << "  PMT count        : " << _pmt_v.size() << std::endl << std::endl << std::flush;
  }


}

#if USING_LARSOFT == 0
namespace flashmatch{

  DetectorSpecs::DetectorSpecs(const Config_t& p) 
    : LoggerFeature("DetectorSpecs")
  {

    auto max_pt = p.get<std::vector<double> >("ActiveVolumeMax");
    auto min_pt = p.get<std::vector<double> >("ActiveVolumeMin");
    auto photon_max_pt = p.get<std::vector<double> >("PhotonLibraryVolumeMax");
    auto photon_min_pt = p.get<std::vector<double> >("PhotonLibraryVolumeMin");
    auto nvoxels = p.get<std::vector<int> >("PhotonLibraryNvoxels");
    auto nopdetchannels = p.get<int>("PhotonLibraryNOpDetChannels");
    assert(max_pt.size() == 3);
    assert(min_pt.size() == 3);
    assert(max_pt[0] >= min_pt[0] &&
     max_pt[1] >= min_pt[1] &&
     max_pt[2] >= min_pt[2]);
    _bbox = geoalgo::AABox(min_pt[0],min_pt[1],min_pt[2],max_pt[0],max_pt[1],max_pt[2]);
    //std::cout<<_bbox.Min()[0]<<" "<<_bbox.Min()[1]<<" "<<_bbox.Min()[2]<<std::endl;
    //std::cout<<_bbox.Max()[0]<<" "<<_bbox.Max()[1]<<" "<<_bbox.Max()[2]<<std::endl;

    assert(photon_max_pt.size() == 3);
    assert(photon_min_pt.size() == 3);
    assert(photon_max_pt[0] >= photon_min_pt[0] &&
           photon_max_pt[1] >= photon_min_pt[1] &&
           photon_max_pt[2] >= photon_min_pt[2]);
    _photon_bbox = geoalgo::AABox(photon_min_pt[0], photon_min_pt[1], photon_min_pt[2], photon_max_pt[0], photon_max_pt[1], photon_max_pt[2]);

    phot::PhotonVisibilityService& photon_library = phot::PhotonVisibilityService::GetME();
    photon_library.LoadLibrary();
    photon_library.SetMaxX(photon_max_pt[0]);
    photon_library.SetMaxY(photon_max_pt[1]);
    photon_library.SetMaxZ(photon_max_pt[2]);
    photon_library.SetMinX(photon_min_pt[0]);
    photon_library.SetMinY(photon_min_pt[1]);
    photon_library.SetMinZ(photon_min_pt[2]);
    photon_library.SetNvoxelsX(nvoxels[0]);
    photon_library.SetNvoxelsY(nvoxels[1]);
    photon_library.SetNvoxelsZ(nvoxels[2]);
    photon_library.SetNOpDetChannels(nopdetchannels);

    size_t ch=0;
    _pmt_v.clear();
    while(1) {
      std::string key = "PMT" + std::to_string(ch);
      if(!p.contains_value(key)) break;
      geoalgo::Point_t pmt(p.get<std::vector<double> >(key));
      assert(pmt.size()==3);
      _pmt_v.push_back(pmt);
      ch++;
    }

    _drift_velocity = p.get<double>("DriftVelocity");
    _light_yield = p.get<double>("LightYield");
    _MIPdEdx = p.get<double>("MIPdEdx");

  }

  const geoalgo::AABox& DetectorSpecs::ActiveVolume(int tpc, int cryo) const {
    auto iter = _bbox_map.find(std::pair<int,int>(tpc, cryo));
    if (iter == _bbox_map.end()) {
      FLASH_CRITICAL() << "Boundary box map doesn't contain cryo " << cryo
                       << " or tpc " << tpc << "!" << std::endl;
      throw OpT0FinderException();
    }
    return iter->second;
  }

  float DetectorSpecs::GetVisibility(double x, double y, double z, unsigned int opch) const
  { return phot::PhotonVisibilityService::GetME().GetVisibility(x,y,z,opch); }

  float DetectorSpecs::GetVisibilityReflected(double x, double y, double z, unsigned int opch) const
  { return -1; }

  float DetectorSpecs::GetVisibility(int vox_id, unsigned int opch) const
  { return phot::PhotonVisibilityService::GetME().GetLibraryEntry(vox_id,opch);}

  float DetectorSpecs::GetVisibilityReflected(int vox_id, unsigned int opch) const
  { return -1; }

  const std::vector<float>& DetectorSpecs::GetLibraryEntries(int vox_id) const
  { return phot::PhotonVisibilityService::GetME().GetLibraryData()[vox_id]; }

  const sim::PhotonVoxelDef& DetectorSpecs::GetVoxelDef() const
  {
    return phot::PhotonVisibilityService::GetME().GetVoxelDef();
  }

}

#else

namespace flashmatch{
  DetectorSpecs::DetectorSpecs(const Config_t& p){
    auto const clock_data = ::art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
    auto const det_prop = ::art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataForJob(clock_data);
    _drift_velocity = det_prop.DriftVelocity();
    if(_drift_velocity != p.get<double>("DriftVelocity",_drift_velocity)) {
      FLASH_CRITICAL() << "Drift velocity is set in the config as " 
      << p.get<double>("DriftVelocity") << " but disagrees with larsoft " << _drift_velocity << std::endl;
      throw OpT0FinderException();
    }
    _light_yield = p.get<double>("LightYield");
    _MIPdEdx = p.get<double>("MIPdEdx");
    std::vector<size_t> cryo_id_v;
    this->EnableCryostats(cryo_id_v);
  }

  void DetectorSpecs::EnableCryostats(std::vector<size_t> cryo_id_v)
  {

    ::art::ServiceHandle<geo::Geometry> const geo;
    _pmt_v.clear();
    _pmt_v.reserve(geo->NOpDets());

    double global_x_min = 1e9, global_x_max = -1e9;
    double global_y_min = 1e9, global_y_max = -1e9;
    double global_z_min = 1e9, global_z_max = -1e9;

    if(cryo_id_v.empty())
      for(size_t i=0; i<geo->Ncryostats(); ++i) cryo_id_v.push_back(i);

    for (auto const& cryo_id : cryo_id_v) {

      auto const& cryo_geo = geo->Cryostat(cryo_id);

      _pmt_v.reserve(_pmt_v.size()+cryo_geo.NOpDet());
    
      for (size_t opdet = 0; opdet < cryo_geo.NOpDet(); opdet++) {

        std::vector<double> pos(3, 0.);
        cryo_geo.OpDet(opdet).GetCenter(&pos[0]);

        geoalgo::Point_t pmt(pos);
        _pmt_v.push_back(pmt);
      }


      for (size_t tpc_id = 0; tpc_id < cryo_geo.NTPC(); tpc_id++) {
        auto const& tpc_geo = cryo_geo.TPC(tpc_id);
        double x_min = tpc_geo.GetCenter().X() - tpc_geo.HalfWidth();
        double x_max = tpc_geo.GetCenter().X() + tpc_geo.HalfWidth();

        double y_min = tpc_geo.GetCenter().Y() - tpc_geo.HalfHeight();
        double y_max = tpc_geo.GetCenter().Y() + tpc_geo.HalfHeight();

        double z_min = tpc_geo.GetCenter().Z() - tpc_geo.HalfLength();
        double z_max = tpc_geo.GetCenter().Z() + tpc_geo.HalfLength();

        if (x_min < global_x_min) global_x_min = x_min;
        if (x_max > global_x_max) global_x_max = x_max;
        if (y_min < global_y_min) global_y_min = y_min;
        if (y_max > global_y_max) global_y_max = y_max;
        if (z_min < global_z_min) global_z_min = z_min;
        if (z_max > global_z_max) global_z_max = z_max;

        auto pair = std::pair<int,int>(tpc_id, cryo_id);
        _bbox_map[pair] = geoalgo::AABox(x_min, y_min, z_min, x_max, y_max, z_max);
      }
    }

    _bbox = geoalgo::AABox(global_x_min, global_y_min, global_z_min,
                           global_x_max, global_y_max, global_z_max);
    
    art::ServiceHandle<phot::PhotonVisibilityService const> pvs;
    auto lower_pt = this->GetVoxelDef().GetRegionLowerCorner<geo::Point_t>();
    auto upper_pt = this->GetVoxelDef().GetRegionUpperCorner<geo::Point_t>();
    _photon_bbox = geoalgo::AABox(lower_pt.X(),lower_pt.Y(),lower_pt.Z(),upper_pt.X(),upper_pt.Y(),upper_pt.Z());
  }

  float DetectorSpecs::GetVisibility(double x, double y, double z, unsigned int opch) const {
    art::ServiceHandle<phot::PhotonVisibilityService const> pvs;
    geo::Point_t pt(x,y,z);
    return pvs->GetVisibility(pt,opch,false);
  }

  float DetectorSpecs::GetVisibilityReflected(double x, double y, double z, unsigned int opch) const {
    art::ServiceHandle<phot::PhotonVisibilityService const> pvs;
    geo::Point_t pt(x,y,z);
    return pvs->GetVisibility(pt,opch,true);
  }

  float DetectorSpecs::GetVisibility(int vox_id, unsigned int opch) const {
    art::ServiceHandle<phot::PhotonVisibilityService const> pvs;
    return pvs->GetLibraryEntry(vox_id,opch,false);
  }

  float DetectorSpecs::GetVisibilityReflected(int vox_id, unsigned int opch) const {
    art::ServiceHandle<phot::PhotonVisibilityService const> pvs;
    return pvs->GetLibraryEntry(vox_id,opch,true);
  }

  const sim::PhotonVoxelDef& DetectorSpecs::GetVoxelDef() const {
    art::ServiceHandle<phot::PhotonVisibilityService const> pvs;
    return pvs->GetVoxelDef();
  }

  phot::IPhotonLibrary::Counts_t DetectorSpecs::GetLibraryEntries(int vox_id, bool reflWanted) const {
    art::ServiceHandle<phot::PhotonVisibilityService const> pvs;
    return pvs->GetLibraryEntries(vox_id,reflWanted);
  }

  const geoalgo::AABox& DetectorSpecs::ActiveVolume(int tpc, int cryo) const {
    auto iter = _bbox_map.find(std::pair<int,int>(tpc, cryo));
    if (iter == _bbox_map.end()) {
      FLASH_CRITICAL() << "Boundary box map doesn't contain cryo " << cryo
                       << " or tpc " << tpc << "!" << std::endl;
      throw OpT0FinderException();
    }
    return iter->second;
  }

}
#endif

#endif

