#include "FiducialVolumeCosmicIdAlg.h"

namespace ana{

FiducialVolumeCosmicIdAlg::FiducialVolumeCosmicIdAlg(const core::ProviderManager &manager, const Config& config){
  this->reconfigure(manager, config);
}


FiducialVolumeCosmicIdAlg::FiducialVolumeCosmicIdAlg(){}


FiducialVolumeCosmicIdAlg::~FiducialVolumeCosmicIdAlg(){}

void FiducialVolumeCosmicIdAlg::reconfigure(const core::ProviderManager &manager, const Config& config){
  fFiducialVolumes.clear();

  fMinX = config.FiducialCuts().MinX(); 
  fMinY = config.FiducialCuts().MinY(); 
  fMinZ = config.FiducialCuts().MinZ(); 
  fMaxX = config.FiducialCuts().MaxX(); 
  fMaxY = config.FiducialCuts().MaxY(); 
  fMaxZ = config.FiducialCuts().MaxZ();

  // get the set of fiducial volumes
  for (auto const &cryo: manager.GetGeometryProvider()->IterateCryostats()) {
    geo::GeometryCore::TPC_iterator tstart = manager.GetGeometryProvider()->begin_TPC(cryo.ID()),
                                    tend = manager.GetGeometryProvider()->end_TPC(cryo.ID());

    double Min_X = std::min_element(tstart, tend, [](auto &lhs, auto &rhs) { return lhs.MinX() < rhs.MinX(); })->MinX();
    double Min_Y = std::min_element(tstart, tend, [](auto &lhs, auto &rhs) { return lhs.MinY() < rhs.MinY(); })->MinY();
    double Min_Z = std::min_element(tstart, tend, [](auto &lhs, auto &rhs) { return lhs.MinZ() < rhs.MinZ(); })->MinZ();
    
    double Max_X = std::max_element(tstart, tend, [](auto &lhs, auto &rhs) { return lhs.MaxX() < rhs.MaxX(); })->MaxX();
    double Max_Y = std::max_element(tstart, tend, [](auto &lhs, auto &rhs) { return lhs.MaxY() < rhs.MaxY(); })->MaxY();
    double Max_Z = std::max_element(tstart, tend, [](auto &lhs, auto &rhs) { return lhs.MaxZ() < rhs.MaxZ(); })->MaxZ();

    fFiducialVolumes.emplace_back(Min_X + fMinX, Max_X - fMaxX, Min_Y + fMinY, Max_Y - fMaxY, Min_Z + fMinZ, Max_Z - fMaxZ);
  }

  return;
}

// Check if point in fiducial volume used by this algorithm
bool FiducialVolumeCosmicIdAlg::InFiducial(geo::Point_t point){
  for (const geo::BoxBoundedGeo vol: fFiducialVolumes) {
    if (vol.ContainsPosition(point)) return true;
  }
  return false;
}

// Check both start and end points of track are in fiducial volume
bool FiducialVolumeCosmicIdAlg::FiducialVolumeCosmicId(recob::Track track){
  
  bool startInFiducial = InFiducial(track.Vertex());

  bool endInFiducial = InFiducial(track.End());

  if(!startInFiducial && !endInFiducial)  return true;
  
  return false;
  
}


}
