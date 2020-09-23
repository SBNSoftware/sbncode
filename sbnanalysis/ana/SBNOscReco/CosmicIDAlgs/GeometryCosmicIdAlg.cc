#include "GeometryCosmicIdAlg.h"

namespace ana{

GeometryCosmicIdAlg::GeometryCosmicIdAlg(const core::ProviderManager &manager, const Config& config){
  fGeometry = manager.GetGeometryProvider();

  this->reconfigure(config);

}


GeometryCosmicIdAlg::GeometryCosmicIdAlg(){

}


GeometryCosmicIdAlg::~GeometryCosmicIdAlg(){

}


void GeometryCosmicIdAlg::reconfigure(const Config& config){

  return;
}

// Remove any tracks in different TPC to beam activity
bool GeometryCosmicIdAlg::GeometryCosmicId(recob::Track &track, std::vector<art::Ptr<recob::Hit>> &hits, std::map<geo::TPCID, bool> &tpc_flashes) {

  // Remove any tracks that are detected in one TPC and reconstructed in another
  geo::TPCID tpcid =  CosmicIdUtils::DetectedInTPC(hits); 
  if (!tpcid) return true;

  geo::TPCGeo tpc = fGeometry->GetElement(tpcid);

  geo::Point_t start = track.Start();
  geo::Point_t end = track.End();

  // Check the start/end points are in same TPC (track shifted into other TPC because time outside of beam)
  if (!tpc.ContainsPosition(start) || !tpc.ContainsPosition(end)) return true;

  // if there was a flash in time with this track, keep it
  return tpc_flashes.at(tpcid);
}


}
