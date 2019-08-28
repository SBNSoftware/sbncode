#include "CrtHitCosmicIdAlg.h"

namespace ana{

CrtHitCosmicIdAlg::CrtHitCosmicIdAlg(const Config& config){

  this->reconfigure(config);

} //CrtHitCosmicIdAlg()


CrtHitCosmicIdAlg::CrtHitCosmicIdAlg(){

} //CrtHitCosmicIdAlg()


CrtHitCosmicIdAlg::~CrtHitCosmicIdAlg(){

} //~CrtHitCosmicIdAlg()


void CrtHitCosmicIdAlg::reconfigure(const Config& config){

  t0Alg = config.T0Alg();
  fBeamTimeMin = config.BeamTimeLimits().BeamTimeMin();
  fBeamTimeMax = config.BeamTimeLimits().BeamTimeMax();

  return;
} //reconfigure()


// Returns true if matched to CRTHit outside beam time
bool CrtHitCosmicIdAlg::CrtHitCosmicId(recob::Track track, std::vector<sbnd::crt::CRTHit> crtHits, const art::Event& event){

  // Get the closest matched time from CRT hits
  double crtHitTime = t0Alg.T0FromCRTHits(track, crtHits, event);

  // If time is valid and outside the beam time then tag as a cosmic
  if(crtHitTime != -99999 && (crtHitTime < fBeamTimeMin || crtHitTime > fBeamTimeMax)) return true;

  return false;

} //CrtHitCosmicId()
 
}
