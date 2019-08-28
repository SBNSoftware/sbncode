#include "CosmicIdAlg.h"

namespace ana{

CosmicIdAlg::CosmicIdAlg(const Config& config){

  this->reconfigure(config);

}


CosmicIdAlg::CosmicIdAlg(){

}


CosmicIdAlg::~CosmicIdAlg(){

}


void CosmicIdAlg::reconfigure(const Config& config){

  fTpcTrackModuleLabel = config.TpcTrackModuleLabel();
  fPandoraLabel        = config.PandoraLabel();
  fCrtHitModuleLabel   = config.CrtHitModuleLabel();
  fCrtTrackModuleLabel = config.CrtTrackModuleLabel();
  fCaloModuleLabel     = config.CaloModuleLabel();

  fApplyFiducialCut  = config.ApplyFiducialCut();
  fApplyStoppingCut  = config.ApplyStoppingCut();
  fApplyGeometryCut  = config.ApplyGeometryCut();
  fApplyCpaCrossCut  = config.ApplyCpaCrossCut();
  fApplyApaCrossCut  = config.ApplyApaCrossCut();
  fApplyCrtTrackCut  = config.ApplyCrtTrackCut();
  fApplyCrtHitCut    = config.ApplyCrtHitCut();
  fApplyPandoraT0Cut = config.ApplyPandoraT0Cut();

  fOriginalSettings.push_back(fApplyFiducialCut);
  fOriginalSettings.push_back(fApplyStoppingCut);
  fOriginalSettings.push_back(fApplyGeometryCut);
  fOriginalSettings.push_back(fApplyCpaCrossCut);
  fOriginalSettings.push_back(fApplyApaCrossCut);
  fOriginalSettings.push_back(fApplyCrtTrackCut);
  fOriginalSettings.push_back(fApplyCrtHitCut);
  fOriginalSettings.push_back(fApplyPandoraT0Cut);

  fUseTrackAngleVeto    = config.UseTrackAngleVeto();
  fMinSecondTrackLength = config.MinSecondTrackLength();
  fMinVertexDistance    = config.MinVertexDistance();
  fMinMergeAngle        = config.MinMergeAngle();

  fvTag = config.FVTagAlg();
  spTag = config.SPTagAlg();
  ccTag = config.CCTagAlg();
  acTag = config.ACTagAlg();
  chTag = config.CHTagAlg();
  ctTag = config.CTTagAlg();
  ptTag = config.PTTagAlg();

  fBeamTimeMin = config.BeamTimeLimits().BeamTimeMin();
  fBeamTimeMax = config.BeamTimeLimits().BeamTimeMax();

  return;
}

// Change which cuts are run
void CosmicIdAlg::SetCuts(bool FV, bool SP, bool Geo, bool CC, bool AC, bool CT, bool CH, bool PT){

  fApplyFiducialCut = FV;
  fApplyStoppingCut = SP;
  fApplyGeometryCut = Geo;
  fApplyCpaCrossCut = CC;
  fApplyApaCrossCut = AC;
  fApplyCrtTrackCut = CT;
  fApplyCrtHitCut = CH;
  fApplyPandoraT0Cut = PT;

}

// Reset which cuts are run from fhicl parameters
void CosmicIdAlg::ResetCuts(){

  fApplyFiducialCut = fOriginalSettings[0];
  fApplyStoppingCut = fOriginalSettings[1];
  fApplyGeometryCut = fOriginalSettings[2];
  fApplyCpaCrossCut = fOriginalSettings[3];
  fApplyApaCrossCut = fOriginalSettings[4];
  fApplyCrtTrackCut = fOriginalSettings[5];
  fApplyCrtHitCut = fOriginalSettings[6];
  fApplyPandoraT0Cut = fOriginalSettings[7];

}

// Run cuts to decide if track looks like a cosmic
bool CosmicIdAlg::CosmicId(recob::Track track, const art::Event& event, std::vector<double> t0Tpc0, std::vector<double> t0Tpc1){

  // Get associations between tracks and hit collections
  auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTpcTrackModuleLabel);
  art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTpcTrackModuleLabel);
  std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(track.ID());
  // Get associations between track and calorimetry collections
  art::FindManyP<anab::Calorimetry> findManyCalo(tpcTrackHandle, event, fCaloModuleLabel);

  // Tag cosmics from pandora T0 associations
  if(fApplyPandoraT0Cut){
    if(ptTag.PandoraT0CosmicId(track, event)) return true;
  }    

  // Tag cosmics which enter and exit the TPC
  if(fApplyFiducialCut){
    if(fvTag.FiducialVolumeCosmicId(track)) return true;
  }

  // Tag cosmics which enter the TPC and stop
  if(fApplyStoppingCut){
    std::vector<art::Ptr<anab::Calorimetry>> calos = findManyCalo.at(track.ID());
    if(spTag.StoppingParticleCosmicId(track, calos)) return true;
  }

  // Tag cosmics in other TPC to beam activity
  if(fApplyGeometryCut){
    bool tpc0Flash = CosmicIdUtils::BeamFlash(t0Tpc0, fBeamTimeMin, fBeamTimeMax);
    bool tpc1Flash = CosmicIdUtils::BeamFlash(t0Tpc1, fBeamTimeMin, fBeamTimeMax);
    if(geoTag.GeometryCosmicId(track, hits, tpc0Flash, tpc1Flash)) return true;
  }

  // Tag cosmics which cross the CPA
  if(fApplyCpaCrossCut){
    std::vector<recob::Track> tracks;
    for(auto const& tpcTrack : (*tpcTrackHandle)){
      tracks.push_back(tpcTrack);
    }

    if(ccTag.CpaCrossCosmicId(track, tracks, findManyHits)) return true;
  }

  // Tag cosmics which cross the APA
  if(fApplyApaCrossCut){
    if(acTag.ApaCrossCosmicId(track, hits, t0Tpc0, t0Tpc1)) return true;
  }

  // Tag cosmics which match CRT tracks
  if(fApplyCrtTrackCut){
    auto crtTrackHandle = event.getValidHandle<std::vector<crt::CRTTrack>>(fCrtTrackModuleLabel);
    std::vector<crt::CRTTrack> crtTracks = (*crtTrackHandle);

    if(ctTag.CrtTrackCosmicId(track, crtTracks, event)) return true;
  }

  // Tag cosmics which match CRT hits
  if(fApplyCrtHitCut){
    auto crtHitHandle = event.getValidHandle<std::vector<crt::CRTHit>>(fCrtHitModuleLabel);
    std::vector<crt::CRTHit> crtHits = (*crtHitHandle);

    if(chTag.CrtHitCosmicId(track, crtHits, event)) return true;
  }

  return false;

}

// Run cuts to decide if PFParticle looks like a cosmic
bool CosmicIdAlg::CosmicId(recob::PFParticle pfparticle, std::map< size_t, art::Ptr<recob::PFParticle> > pfParticleMap, const art::Event& event, std::vector<double> t0Tpc0, std::vector<double> t0Tpc1){

  // Get associations between pfparticles and tracks
  art::Handle< std::vector<recob::PFParticle> > pfParticleHandle;
  event.getByLabel(fPandoraLabel, pfParticleHandle);
  art::FindManyP< recob::Track > pfPartToTrackAssoc(pfParticleHandle, event, fTpcTrackModuleLabel);

  // Get associations between tracks and hits/calorimetry collections
  auto tpcTrackHandle = event.getValidHandle<std::vector<recob::Track>>(fTpcTrackModuleLabel);
  art::FindManyP<recob::Hit> findManyHits(tpcTrackHandle, event, fTpcTrackModuleLabel);
  art::FindManyP<anab::Calorimetry> findManyCalo(tpcTrackHandle, event, fCaloModuleLabel);

  // Loop over all the daughters of the PFParticles and get associated tracks
  std::vector<recob::Track> nuTracks;
  for (const size_t daughterId : pfparticle.Daughters()){
  
    // Get tracks associated with daughter
    art::Ptr<recob::PFParticle> pParticle = pfParticleMap.at(daughterId);
    const std::vector< art::Ptr<recob::Track> > associatedTracks(pfPartToTrackAssoc.at(pParticle.key()));
    if(associatedTracks.size() != 1) continue;

    recob::Track track = *associatedTracks.front();
    nuTracks.push_back(track);
    
    
  }

  // Tag cosmics from pandora T0 associations
  if(fApplyPandoraT0Cut){
    if(ptTag.PandoraT0CosmicId(pfparticle, pfParticleMap, event)) return true;
  }

  // Not a cosmic if there are only showers assiciated with PFParticle
  if(nuTracks.size() == 0) return false;

  // Sort all daughter tracks by length
  std::sort(nuTracks.begin(), nuTracks.end(), [](auto& left, auto& right){
              return left.Length() > right.Length();});

  // Select longest track as the cosmic candidate
  recob::Track track = nuTracks[0];
  std::vector<art::Ptr<recob::Hit>> hits = findManyHits.at(track.ID());

  // Tag cosmics which enter and exit the TPC
  if(fApplyFiducialCut){
    if(fvTag.FiducialVolumeCosmicId(track)) return true;
  }

  // Tag cosmics in other TPC to beam activity
  if(fApplyGeometryCut){
    bool tpc0Flash = CosmicIdUtils::BeamFlash(t0Tpc0, fBeamTimeMin, fBeamTimeMax);
    bool tpc1Flash = CosmicIdUtils::BeamFlash(t0Tpc1, fBeamTimeMin, fBeamTimeMax);
    if(geoTag.GeometryCosmicId(track, hits, tpc0Flash, tpc1Flash)) return true;
  }

  // Tag cosmics which match CRT tracks
  if(fApplyCrtTrackCut){
    auto crtTrackHandle = event.getValidHandle<std::vector<crt::CRTTrack>>(fCrtTrackModuleLabel);
    std::vector<crt::CRTTrack> crtTracks = (*crtTrackHandle);

    if(ctTag.CrtTrackCosmicId(track, crtTracks, event)) return true;
  }

  // Tag cosmics which cross the CPA
  if(fApplyCpaCrossCut){
    std::vector<recob::Track> tracks;
    for(auto const& tpcTrack : (*tpcTrackHandle)){
      tracks.push_back(tpcTrack);
    }

    if(ccTag.CpaCrossCosmicId(track, tracks, findManyHits)) return true;
  }

  // Find second longest particle if trying to merge tracks
  std::vector<std::pair<recob::Track, double>> secondaryTracks;
  if(fUseTrackAngleVeto && nuTracks.size() > 1){
    TVector3 start = track.Vertex<TVector3>();
    TVector3 end = track.End<TVector3>();

    // Loop over the secondary tracks
    // Find smallest angle between primary track and any secondary tracks above a certain length
    for(size_t i = 1; i < nuTracks.size(); i++){
      recob::Track track2 = nuTracks[i];
      // Only consider secondary tracks longer than some limit (try to exclude michel electrons)
      if(track2.Length() < fMinSecondTrackLength) continue;
      TVector3 start2 = track2.Vertex<TVector3>();
      TVector3 end2 = track2.End<TVector3>();
      // Do they share the same vertex? (no delta rays)
      if((start-start2).Mag() < fMinVertexDistance){ 
        double angle = (end - start).Angle(end2 - start2);
        secondaryTracks.push_back(std::make_pair(track2, angle));
      }
    }
  }

  // If there is a valid secondary track
  if(secondaryTracks.size() > 0){
    // Sort tracks by smallest angle
    std::sort(secondaryTracks.begin(), secondaryTracks.end(), [](auto& left, auto& right){
              return left.second < right.second;});
    // If secondary track angle is compatible with split track (near 180) then try to merge
    if(secondaryTracks[0].second > fMinMergeAngle){
      recob::Track track2 = secondaryTracks[0].first;

      // Check fiducial volume containment assuming merged track
      if(fApplyFiducialCut){
        if(!fvTag.InFiducial(track.End()) && !fvTag.InFiducial(track2.End())) return true;
      }

      // Check if stopping applies to merged track
      if(fApplyStoppingCut){
        // Apply stopping cut to the longest track
        std::vector<art::Ptr<anab::Calorimetry>> calos = findManyCalo.at(track.ID());
        if(spTag.StoppingParticleCosmicId(track, calos)) return true;
        // Apply stopping cut assuming the tracks are split
        std::vector<art::Ptr<anab::Calorimetry>> calos2 = findManyCalo.at(track2.ID());
        if(spTag.StoppingParticleCosmicId(track, track2, calos, calos2)) return true;
      }

      // Check if either track crosses APA
      if(fApplyApaCrossCut){
        // Apply apa crossing cut to the longest track
        if(acTag.ApaCrossCosmicId(track, hits, t0Tpc0, t0Tpc1)) return true;
        // Also apply to secondary track FIXME need to check primary track doesn't go out of bounds
        std::vector<art::Ptr<recob::Hit>> hits2 = findManyHits.at(track2.ID());
        if(acTag.ApaCrossCosmicId(track2, hits2, t0Tpc0, t0Tpc1)) return true;
      }

      // Check if either track matches CRT hit
      if(fApplyCrtHitCut){
        // Apply crt hit match cut to both tracks
        auto crtHitHandle = event.getValidHandle<std::vector<crt::CRTHit>>(fCrtHitModuleLabel);
        std::vector<crt::CRTHit> crtHits = (*crtHitHandle);
        if(chTag.CrtHitCosmicId(track, crtHits, event)) return true;
        if(chTag.CrtHitCosmicId(track2, crtHits, event)) return true;
      }
    }
    // Don't apply other cuts if angle between tracks is consistent with neutrino interaction
  }

  // If there's only one track from vertex then apply cuts as normal
  if(secondaryTracks.size() == 0){

    // Tag cosmics which enter the TPC and stop
    if(fApplyStoppingCut){
      std::vector<art::Ptr<anab::Calorimetry>> calos = findManyCalo.at(track.ID());
      if(spTag.StoppingParticleCosmicId(track, calos)) return true;
    }

    // Tag cosmics which cross the APA
    if(fApplyApaCrossCut){
      if(acTag.ApaCrossCosmicId(track, hits, t0Tpc0, t0Tpc1)) return true;
    }

    // Tag cosmics which match CRT hits
    if(fApplyCrtHitCut){
      auto crtHitHandle = event.getValidHandle<std::vector<crt::CRTHit>>(fCrtHitModuleLabel);
      std::vector<crt::CRTHit> crtHits = (*crtHitHandle);

      if(chTag.CrtHitCosmicId(track, crtHits, event)) return true;
    }
  }

  return false;

}


}
