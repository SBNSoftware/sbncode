#include "CrtTrackCosmicIdAlg.h"

namespace ana{

CrtTrackCosmicIdAlg::CrtTrackCosmicIdAlg(const core::ProviderManager &manager, const Config& config){

  this->reconfigure(manager, config);

}


CrtTrackCosmicIdAlg::CrtTrackCosmicIdAlg(){

}


CrtTrackCosmicIdAlg::~CrtTrackCosmicIdAlg(){

}


void CrtTrackCosmicIdAlg::reconfigure(const core::ProviderManager &manager, const Config& config){

  trackMatchAlg = sbnd::CRTTrackMatchAlg(config.TrackMatchAlg(), manager.GetGeometryProvider(), (const detinfo::DetectorProperties*)manager.GetDetectorPropertiesProvider());
  fBeamTimeMin = config.BeamTimeLimits().BeamTimeMin();
  fBeamTimeMax = config.BeamTimeLimits().BeamTimeMax();

  return;
}


// Tags track as cosmic if it matches a CRTTrack
bool CrtTrackCosmicIdAlg::CrtTrackCosmicId(recob::Track track, std::vector<art::Ptr<recob::Hit>> hits, std::vector<sbnd::crt::CRTTrack> crtTracks) {

  // Get the closest matching CRT track ID
  int crtID = trackMatchAlg.GetMatchedCRTTrackId(track, hits, crtTracks);

  // If matching failed
  if(crtID == -99999) return false;

  // If track matched to a through going CRT track then it is a cosmic
  if(crtTracks.at(crtID).complete) return true;

  // If it matches a track through just the top planes make sure it is outside of the beam time
  double crtTime = ((double)(int)crtTracks.at(crtID).ts1_ns) * 1e-3; // [us]
  if(crtTime < fBeamTimeMin || crtTime > fBeamTimeMax) return true;

  return false;

}
 
}
