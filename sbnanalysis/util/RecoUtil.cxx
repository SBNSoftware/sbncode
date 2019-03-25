#include "RecoUtil.hh"

namespace util {
  
  int TrueParticleID(const art::Ptr<recob::Hit> hit, core::ProviderManager* provider_manager, bool rollup_unsaved_ids){
    std::map<int,double> id_to_energy_map;
    std::vector<sim::TrackIDE> track_ides = provider_manager->GetBackTrackerProvider()->HitToTrackIDEs(hit);
    for (unsigned int idIt = 0; idIt < track_ides.size(); ++idIt) {
      int id = track_ides.at(idIt).trackID;
      if (rollup_unsaved_ids) id = std::abs(id);
      double energy = track_ides.at(idIt).energy;
      id_to_energy_map[id]+=energy;
    }
    //Now loop over the map to find the maximum contributor
    double likely_particle_contrib_energy = -99999;
    int likely_track_id = 0;
    for (std::map<int,double>::iterator mapIt = id_to_energy_map.begin(); 
         mapIt != id_to_energy_map.end(); mapIt++){ 
      double particle_contrib_energy = mapIt->second;
      if (particle_contrib_energy > likely_particle_contrib_energy){
        likely_particle_contrib_energy = particle_contrib_energy;
        likely_track_id = mapIt->first;
      }
    }
    return likely_track_id;
  } // TrueParticleID
   
  int TrueParticleIDFromTotalRecoHits(const std::vector< art::Ptr< recob::Hit > >& hits, core::ProviderManager* provider_manager, bool rollup_unsaved_ids) {
    // Make a map of the tracks which are associated with this object and the number of hits they are the primary contributor to
    std::map< int, int > trackMap;
    for (std::vector< art::Ptr< recob::Hit > >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
      art::Ptr< recob::Hit > hit = *hitIt;
      int trackID = TrueParticleID(hit, provider_manager, rollup_unsaved_ids);
      trackMap[trackID]++;
    }

    // Pick the track which is the primary contributor to the most hits as the 'true track'
    int objectTrack = -99999;
    int highestCount = -1;
    int NHighestCounts = 0;
    for (std::map< int, int >::iterator trackIt = trackMap.begin(); trackIt != trackMap.end(); ++trackIt) {
      if (trackIt->second > highestCount) {
        highestCount = trackIt->second;
        objectTrack  = trackIt->first;
        NHighestCounts = 1;
      }
      else if (trackIt->second == highestCount){
        NHighestCounts++;
      }
    }
    if (NHighestCounts > 1){
      std::cout << "util::TrueParticleIDFromTotalRecoHits - There are " << NHighestCounts << " particles which tie for highest number of contributing hits ( " << highestCount << " hits). Using util::TrueParticleIDFromTotalTrueEnergy instead." << std::endl;
      objectTrack = TrueParticleIDFromTotalTrueEnergy(hits, provider_manager, rollup_unsaved_ids);
    }
    return objectTrack;
  } // TrueParticleIDFromTotalRecoHits

  int TrueParticleIDFromTotalTrueEnergy(const std::vector< art::Ptr< recob::Hit > >& hits, core::ProviderManager* provider_manager, bool rollup_unsaved_ids) {
    std::map<int,double> trackIDToEDepMap;
    for (std::vector< art::Ptr< recob::Hit > >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
      art::Ptr< recob::Hit > hit = *hitIt;
      std::vector< sim::TrackIDE > trackIDs = provider_manager->GetBackTrackerProvider()->HitToTrackIDEs(hit);
      for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
        int id = trackIDs[idIt].trackID; 
        if (rollup_unsaved_ids) id = std::abs(id);
        trackIDToEDepMap[id] += trackIDs[idIt].energy;
      } 
    } 

    // Loop over the map and find the track which contributes the highest energy to the hit vector
    double maxenergy = -1;
    int objectTrack = -99999;
    for (std::map< int, double >::iterator mapIt = trackIDToEDepMap.begin(); mapIt != trackIDToEDepMap.end(); mapIt++){
      double energy = mapIt->second;
      double trackid = mapIt->first;
      if (energy > maxenergy){
        maxenergy = energy;
        objectTrack = trackid;
      }
    }
    return objectTrack;
  } // TrueParticleFromTotalTrueEnergy
} // namespace util

