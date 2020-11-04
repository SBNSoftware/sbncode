#include "CrtHitCosmicIdAlg.h"

namespace ana{

CrtHitCosmicIdAlg::CrtHitCosmicIdAlg(const core::ProviderManager &manager, const Config& config){

  this->reconfigure(manager, config);

} //CrtHitCosmicIdAlg()


CrtHitCosmicIdAlg::CrtHitCosmicIdAlg(){

} //CrtHitCosmicIdAlg()


CrtHitCosmicIdAlg::~CrtHitCosmicIdAlg(){

} //~CrtHitCosmicIdAlg()


void CrtHitCosmicIdAlg::reconfigure(const core::ProviderManager &manager, const Config& config){

  t0Alg = sbnd::CRTT0MatchAlg(config.T0Alg(), manager.GetGeometryProvider(), (const detinfo::DetectorProperties*)manager.GetDetectorPropertiesProvider());
  fBeamTimeMin = config.BeamTimeLimits().BeamTimeMin();
  fBeamTimeMax = config.BeamTimeLimits().BeamTimeMax();

  return;
} //reconfigure()


// Returns true if matched to CRTHit outside beam time
bool CrtHitCosmicIdAlg::CrtHitCosmicId(recob::Track track, std::vector<art::Ptr<recob::Hit>> hits, std::vector<sbn::crt::CRTHit> crtHits) {

  // Get the closest matched time from CRT hits
  double crtHitTime = t0Alg.T0FromCRTHits(track, hits, crtHits);

  // If time is valid and outside the beam time then tag as a cosmic
  if(crtHitTime != -99999 && (crtHitTime < fBeamTimeMin || crtHitTime > fBeamTimeMax)) return true;

  return false;

} //CrtHitCosmicId()
 
}
