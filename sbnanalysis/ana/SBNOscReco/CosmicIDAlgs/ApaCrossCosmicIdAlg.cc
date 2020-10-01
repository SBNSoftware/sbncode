#include "ApaCrossCosmicIdAlg.h"

namespace ana {

ApaCrossCosmicIdAlg::ApaCrossCosmicIdAlg(const core::ProviderManager &manager, const Config& config){
  this->reconfigure(manager, config);
}


ApaCrossCosmicIdAlg::ApaCrossCosmicIdAlg() {
  fDetectorProperties = NULL;
}

ApaCrossCosmicIdAlg::~ApaCrossCosmicIdAlg(){}

void ApaCrossCosmicIdAlg::reconfigure(const core::ProviderManager &manager, const Config& config){
  fDetectorProperties = manager.GetDetectorPropertiesProvider();
  fGeometry = manager.GetGeometryProvider();
  fDistanceLimit = config.DistanceLimit(); 
  fMaxApaDistance = config.MaxApaDistance();
  fBeamTimeMin = config.BeamTimeLimits().BeamTimeMin();
  fBeamTimeMax = config.BeamTimeLimits().BeamTimeMax();

  return;
}


// Get the minimum distance from track to APA for different times
std::pair<double, double> ApaCrossCosmicIdAlg::MinApaDistance(const recob::Track &track, std::vector<double> &t0List, geo::TPCID &tpcid){
  double crossTime = -99999;

  double minDist = 99999;
  double startX = track.Vertex().X();
  double endX = track.End().X();
  geo::Point_t point = track.Vertex();


  // if TPCID is invalid, return null values
  if (!tpcid) {
    return std::make_pair(minDist, crossTime);
  }

  // get the TPC drift direction
  const geo::TPCGeo &tpc = fGeometry->GetElement(tpcid);
  int drift_direction_x = tpc.DriftDirection() == geo::driftdir::kPos ? 1 : -1;

  // get the location of the APA in the TPC -- should be the x-position downstream of the drift direction
  double apa_X = (drift_direction_x == 1) ? tpc.MaxX() : tpc.MinX();

  // Don't try to shift tracks near the APA (may give artificially small distances)
  if (std::abs(startX - apa_X) < fMaxApaDistance || std::abs(endX - apa_X) < fMaxApaDistance)
    return std::make_pair(minDist, crossTime);

  // if drift direction positive use start/end with highest X
  if (drift_direction_x == 1) {
    point = endX < startX ? track.End() : track.Vertex();
  }
  else {
    point = endX > startX ? track.End() : track.Vertex();
  }

  // Shift track by all t0's
  for(auto const& t0 : t0List){
    // If particle crosses the APA before t = 0 the crossing point won't be reconstructed
    if(t0 < 0) continue;
    double shift = t0 * fDetectorProperties->DriftVelocity() * drift_direction_x;
    double shiftedX = point.X() + shift;

    //Check track still in TPC
    if (shiftedX - apa_X > fDistanceLimit) continue;
    // Calculate distance between start/end and APA
    double dist = std::abs(shiftedX - apa_X);
    if(dist < minDist) {
      minDist = dist;
      crossTime = t0;
    }
  }

  return std::make_pair(minDist, crossTime);

}


// Get time by matching tracks which cross the APA
double ApaCrossCosmicIdAlg::T0FromApaCross(const recob::Track &track, std::vector<art::Ptr<recob::Hit>> hits, std::map<geo::CryostatID, std::vector<double>> &t_zeros) {
  // Determine the TPC from hit collection
  geo::TPCID tpc = CosmicIdUtils::DetectedInTPC(hits);
  if (tpc) {
    // Get the minimum distance to the APA and corresponding time
    std::pair<double, double> min = MinApaDistance(track, t_zeros.at(tpc), tpc);
    // Check the distance is within allowed limit
    if(min.first < fDistanceLimit) return min.second;
  }

  return -99999;

}

// Get the distance from track to APA at fixed time
double ApaCrossCosmicIdAlg::ApaDistance(recob::Track track, double t0, std::vector<art::Ptr<recob::Hit>> hits){

  std::vector<double> t0List {t0};
  // Determine the TPC from hit collection
  geo::TPCID tpc = CosmicIdUtils::DetectedInTPC(hits);
  // Get the distance to the APA at the given time
  std::pair<double, double> min = MinApaDistance(track, t0List, tpc);
  return min.first;

}

// Tag tracks with times outside the beam
bool ApaCrossCosmicIdAlg::ApaCrossCosmicId(const recob::Track &track, std::vector<art::Ptr<recob::Hit>> &hits, std::map<geo::CryostatID, std::vector<double>> &t_zeros) {
  double crossTimeThrough = T0FromApaCross(track, hits, t_zeros);
  // If the matched time is outside of the beam time then tag as a cosmic
  if(crossTimeThrough != -99999 && (crossTimeThrough < fBeamTimeMin || crossTimeThrough > fBeamTimeMax)) return true;
  return false;
}

}
