#include "CpaCrossCosmicIdAlg.h"

namespace ana{

CpaCrossCosmicIdAlg::CpaCrossCosmicIdAlg(const core::ProviderManager &manager, const Config& config){
  this->reconfigure(config);
  fDetectorProperties = manager.GetDetectorPropertiesProvider();
}


CpaCrossCosmicIdAlg::CpaCrossCosmicIdAlg(){
  fDetectorProperties = NULL;
}

CpaCrossCosmicIdAlg::~CpaCrossCosmicIdAlg() {}

void CpaCrossCosmicIdAlg::reconfigure(const Config& config){

  fCpaStitchDistance = config.CpaStitchDistance(); 
  fCpaStitchAngle = config.CpaStitchAngle();
  fCpaXDifference = config.CpaXDifference();
  fMinX = config.FiducialCuts().MinX(); 
  fMinY = config.FiducialCuts().MinY(); 
  fMinZ = config.FiducialCuts().MinZ(); 
  fMaxX = config.FiducialCuts().MaxX(); 
  fMaxY = config.FiducialCuts().MaxY(); 
  fMaxZ = config.FiducialCuts().MaxZ();
  fBeamTimeMin = config.BeamTimeLimits().BeamTimeMin();
  fBeamTimeMax = config.BeamTimeLimits().BeamTimeMax();

  return;
}

// Calculate the time by stitching tracks across the CPA
std::pair<double, bool> CpaCrossCosmicIdAlg::T0FromCpaStitching(recob::Track t1, std::vector<recob::Track> tracks){
  
  std::vector<std::pair<double, std::pair<double, bool>>> matchCandidates;
  double matchedTime = -99999;
  std::pair<double, bool> returnVal = std::make_pair(matchedTime, false);

  TVector3 trk1Front = t1.Vertex<TVector3>();
  TVector3 trk1Back = t1.End<TVector3>();
  double closestX1 = std::min(std::abs(trk1Front.X()), std::abs(trk1Back.X()));

  // Loop over all tracks in other TPC
  for(auto & track : tracks){

    TVector3 trk2Front = track.Vertex<TVector3>();
    TVector3 trk2Back = track.End<TVector3>();
    double closestX2 = std::min(std::abs(trk2Front.X()), std::abs(trk2Back.X()));

    // Try to match if their ends have similar x positions
    if(std::abs(closestX1-closestX2) > fCpaXDifference) continue;

    // Find which point is closest to CPA
    TVector3 t1Pos = trk1Front;
    TVector3 t1PosEnd = trk1Back;
    TVector3 t1Dir = t1.VertexDirection<TVector3>();
    if(std::abs(trk1Back.X()) == closestX1){ 
      t1Pos = trk1Back;
      t1PosEnd = trk1Front;
      t1Dir = t1.EndDirection<TVector3>();
    }

    TVector3 t2Pos = trk2Front;
    TVector3 t2PosEnd = trk2Back;
    TVector3 t2Dir = track.VertexDirection<TVector3>();
    if(std::abs(trk2Back.X()) == closestX2){ 
      t2Pos = trk2Back;
      t2PosEnd = trk2Front;
      t2Dir = track.EndDirection<TVector3>();
    }

    // Calculate the angle between the tracks
    double trkCos = std::abs(t1Dir.Dot(t2Dir));
    // Calculate the distance between the tracks at the middle of the CPA
    t1Pos[0] = 0.;
    t2Pos[0] = 0.;
    double dist = (t1Pos-t2Pos).Mag();

    // Does the track enter and exit the fiducial volume when merged?
    geo::Point_t mergeStart {t1PosEnd.X(), t1PosEnd.Y(), t1PosEnd.Z()};
    geo::Point_t mergeEnd {t2PosEnd.X(), t2PosEnd.Y(), t2PosEnd.Z()};
    bool exits = false;
    if(!fTpcGeo.InFiducial(mergeStart, fMinX, fMinY, fMinZ, fMaxX, fMaxY, fMaxZ) 
       && !fTpcGeo.InFiducial(mergeEnd, fMinX, fMinY, fMinZ, fMaxX, fMaxY, fMaxZ)) exits = true;

    // If the distance and angle are within the acceptable limits then record candidate
    if(dist < fCpaStitchDistance && trkCos > cos(TMath::Pi() * fCpaStitchAngle / 180.)){ 
      matchCandidates.push_back(std::make_pair(trkCos, std::make_pair(closestX1, exits)));
    }
  }

  // Choose the candidate with the smallest angle
  if(matchCandidates.size() > 0){
    std::sort(matchCandidates.begin(), matchCandidates.end(), [](auto& left, auto& right){
              return left.first < right.first;});
    double shiftX = matchCandidates[0].second.first;
    matchedTime = -((shiftX - fTpcGeo.CpaWidth())/fDetectorProperties->DriftVelocity()); //subtract CPA width
    returnVal = std::make_pair(matchedTime, matchCandidates[0].second.second);
  }

  return returnVal;
}

// Tag tracks as cosmics from CPA stitching t0
bool CpaCrossCosmicIdAlg::CpaCrossCosmicId(recob::Track track, std::vector<recob::Track> tracks, art::FindManyP<recob::Hit> hitAssoc){

  // Sort tracks by tpc
  std::vector<recob::Track> tpcTracksTPC0;
  std::vector<recob::Track> tpcTracksTPC1;
  // Loop over the tpc tracks
  for(auto const& tpcTrack : tracks){
    // Work out where the associated wire hits were detected
    std::vector<art::Ptr<recob::Hit>> hits = hitAssoc.at(tpcTrack.ID());
    int tpc = fTpcGeo.DetectedInTPC(hits);
    double startX = tpcTrack.Start().X();
    double endX = tpcTrack.End().X();
    if(tpc == 0 && !(startX>0 || endX>0)) tpcTracksTPC0.push_back(tpcTrack);
    else if(tpc == 1 && !(startX<0 || endX<0)) tpcTracksTPC1.push_back(tpcTrack);
  }

  std::vector<art::Ptr<recob::Hit>> hits = hitAssoc.at(track.ID());
  int tpc = fTpcGeo.DetectedInTPC(hits);

  double stitchTime = -99999;
  bool stitchExit = false;
  // Try to match tracks from CPA crossers
  if(tpc == 0){
    std::pair<double, bool> stitchResults = T0FromCpaStitching(track, tpcTracksTPC1);
    stitchTime = stitchResults.first;
    stitchExit = stitchResults.second;
  }
  else if(tpc == 1){
    std::pair<double, bool> stitchResults = T0FromCpaStitching(track, tpcTracksTPC0);
    stitchTime = stitchResults.first;
    stitchExit = stitchResults.second;
  }

  // If tracks are stitched, get time and remove any outside of beam window
  if(stitchTime != -99999 && (stitchTime < fBeamTimeMin || stitchTime > fBeamTimeMax || stitchExit)) return true;
  
  return false;

}

}
