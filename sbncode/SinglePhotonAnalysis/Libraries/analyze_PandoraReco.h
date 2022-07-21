#ifndef SBNCODE_SINGLEPHOTONANALYSIS_ANALYZE_PANDORARECO_H
#define SBNCODE_SINGLEPHOTONANALYSIS_ANALYZE_PANDORARECO_H

#include "sbncode/SinglePhotonAnalysis/Libraries/variables.h"
#include "sbncode/SinglePhotonAnalysis/HelperFunctions/helper_PandoraPFParticles.h"

namespace single_photon
{



  //Analyze Tracks
  void AnalyzeTracks(
      std::vector<PandoraPFParticle> all_PPFPs,
      const std::vector<art::Ptr<recob::Track>>& tracks,
      std::map<art::Ptr<recob::PFParticle>, std::vector<art::Ptr<recob::SpacePoint>>> & pfParticleToSpacePointsMap, 
      std::map<int, art::Ptr<simb::MCParticle> > & MCParticleToTrackIdMap,
      std::map<int, double> &sliceIdToNuScoreMap);

  void AnalyzeTrackCalo(const std::vector<art::Ptr<recob::Track>> &tracks, std::vector<PandoraPFParticle> all_PPFPs);

  void CollectPID(std::vector<art::Ptr<recob::Track>> & tracks,
      std::vector<PandoraPFParticle> all_PPFPs);


  //Analyze falshes
  void AnalyzeFlashes(const std::vector<art::Ptr<recob::OpFlash>>& flashes, art::Handle<std::vector<sbn::crt::CRTHit>> crthit_h, double evt_timeGPS_nsec,  std::map<art::Ptr<recob::OpFlash>, std::vector< art::Ptr<sbn::crt::CRTHit>>> crtvetoToFlashMap);


  //Analyze Showers
  void AnalyzeShowers(
      std::vector<PandoraPFParticle> all_PPFPs,
      const std::vector<art::Ptr<recob::Shower>>& showers,  
      std::map<art::Ptr<recob::Cluster>,  std::vector<art::Ptr<recob::Hit>> >  & clusterToHitMap , 
      double triggeroffset,
      detinfo::DetectorPropertiesData const & theDetector
    );
}
#endif // SBNCODE_SINGLEPHOTONANALYSIS_ANALYZE_PANDORARECO_H
