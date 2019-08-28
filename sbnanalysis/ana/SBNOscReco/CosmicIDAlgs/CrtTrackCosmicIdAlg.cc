#include "CrtTrackCosmicIdAlg.h"

namespace ana{

CrtTrackCosmicIdAlg::CrtTrackCosmicIdAlg(const Config& config){

  this->reconfigure(config);

}


CrtTrackCosmicIdAlg::CrtTrackCosmicIdAlg(){

}


CrtTrackCosmicIdAlg::~CrtTrackCosmicIdAlg(){

}


void CrtTrackCosmicIdAlg::reconfigure(const Config& config){

  trackMatchAlg = config.TrackMatchAlg();
  fBeamTimeMin = config.BeamTimeLimits().BeamTimeMin();
  fBeamTimeMax = config.BeamTimeLimits().BeamTimeMax();

  return;
}


// Tags track as cosmic if it matches a CRTTrack
bool CrtTrackCosmicIdAlg::CrtTrackCosmicId(recob::Track track, std::vector<sbnd::crt::CRTTrack> crtTracks, const art::Event& event){

  // Get the closest matching CRT track ID
  int crtID = trackMatchAlg.GetMatchedCRTTrackId(track, crtTracks, event);

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
