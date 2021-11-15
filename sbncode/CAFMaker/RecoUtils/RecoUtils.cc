#include "RecoUtils.h"

std::vector<std::pair<int, float>> CAFRecoUtils::AllTrueParticleIDEnergyMatches(const detinfo::DetectorClocksData &clockData, const std::vector<art::Ptr<recob::Hit> >& hits, bool rollup_unsaved_ids) {
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;
  std::map<int, float> trackIDToEDepMap;
  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
    art::Ptr<recob::Hit> hit = *hitIt;
    std::vector<sim::TrackIDE> trackIDs = bt_serv->HitToTrackIDEs(clockData, hit);
    for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
      int id = trackIDs[idIt].trackID;
      if (rollup_unsaved_ids) id = std::abs(id);
      id = GetShowerPrimary(id);
      trackIDToEDepMap[id] += trackIDs[idIt].energy;
    }
  }

  std::vector<std::pair<int, float>> ret;
  for (auto const &pair: trackIDToEDepMap) {
    ret.push_back(pair);
  }
  return ret;
}

float CAFRecoUtils::TrackPurity(const detinfo::DetectorClocksData &clockData, int mcparticle_id, const std::vector<art::Ptr<recob::Hit>> &reco_track_hits) {
  // get handle to back tracker
  art::ServiceHandle<cheat::BackTrackerService> bt;

  // get all the hits of the reco track that match the truth track
  const std::vector<art::Ptr<recob::Hit>> matched_reco_track_hits = bt->TrackIdToHits_Ps(clockData, mcparticle_id, reco_track_hits);

  // for each of the hits get the energy coming from the track
  float matched_reco_energy = 0.;
  float reco_energy = 0.;
  for (auto const &matched_reco_track_hit: matched_reco_track_hits) {
    std::vector<sim::IDE> this_hit_IDEs = bt->HitToAvgSimIDEs(clockData, *matched_reco_track_hit);
    for (auto const &ide: this_hit_IDEs) {
      reco_energy += ide.energy;
      if (ide.trackID == mcparticle_id) {
        matched_reco_energy += ide.energy;
      }
    }
  }

  return (reco_energy > 1e-6) ? matched_reco_energy / reco_energy : 1.; 
}

float CAFRecoUtils::TotalHitEnergy(const detinfo::DetectorClocksData &clockData, const std::vector<art::Ptr<recob::Hit> >& hits) {
  art::ServiceHandle<cheat::BackTrackerService> bt_serv;

  float ret = 0.;

  for (std::vector<art::Ptr<recob::Hit> >::const_iterator hitIt = hits.begin(); hitIt != hits.end(); ++hitIt) {
    art::Ptr<recob::Hit> hit = *hitIt;
    std::vector<sim::TrackIDE> trackIDs = bt_serv->HitToTrackIDEs(clockData, hit);
    for (unsigned int idIt = 0; idIt < trackIDs.size(); ++idIt) {
      ret += trackIDs[idIt].energy;
    }
  }
  return ret;
}

float CAFRecoUtils::TrackCompletion(const detinfo::DetectorClocksData &clockData, int mcparticle_id, const std::vector<art::Ptr<recob::Hit>> &reco_track_hits) {
  // get handle to back tracker
  art::ServiceHandle<cheat::BackTrackerService> bt;

  // get all the IDE's of the truth track
  const std::vector<const sim::IDE*> mcparticle_ides = bt->TrackIdToSimIDEs_Ps(mcparticle_id);
  // sum it up
  float mcparticle_energy = 0.;
  for (auto const &ide: mcparticle_ides) {
    mcparticle_energy += ide->energy;
  }

  // get all the hits of the reco track that match the truth track
  const std::vector<art::Ptr<recob::Hit>> matched_reco_track_hits = bt->TrackIdToHits_Ps(clockData, mcparticle_id, reco_track_hits);

  // for each of the hits get the energy coming from the track
  float matched_reco_energy = 0.;
  for (auto const &matched_reco_track_hit: matched_reco_track_hits) {
    std::vector<sim::IDE> this_hit_IDEs = bt->HitToAvgSimIDEs(clockData, *matched_reco_track_hit);
    for (auto const &ide: this_hit_IDEs) {
      if (ide.trackID == mcparticle_id) {
        matched_reco_energy += ide.energy;
      }
    }
  }

  return matched_reco_energy / mcparticle_energy;
}

int CAFRecoUtils::GetShowerPrimary(const int g4ID)
{
  art::ServiceHandle<cheat::ParticleInventoryService> particleInventory;
  const sim::ParticleList& particles = particleInventory->ParticleList();
  const sim::ParticleList::const_iterator part_iter = particles.find(g4ID);
  if(part_iter == particles.end()) return g4ID;

  auto temp_iter = part_iter;
  int primary_id = part_iter->second->TrackId();

  while (std::abs(temp_iter->second->PdgCode()) == 11 || temp_iter->second->PdgCode() == 22)
    {
      primary_id = temp_iter->second->TrackId();
      temp_iter = particles.find(temp_iter->second->Mother());
      if(temp_iter == particles.end()) break;
    }

  return primary_id;
}
