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

std::map<int, std::vector<sbn::ReadoutIDE>> sbn::PrepSimChannels(const std::vector<art::Ptr<sim::SimChannel>> &simchannels, const geo::WireReadoutGeom &wireReadout) {
    std::map<int, std::vector<sbn::ReadoutIDE>> ret;

    for (const art::Ptr<sim::SimChannel> sc : simchannels) {
      // Lookup the wire of this channel
      raw::ChannelID_t channel = sc->Channel();
      std::vector<geo::WireID> maybewire = wireReadout.ChannelToWire(channel);
      geo::WireID thisWire; // Default constructor makes invalid wire
      if (maybewire.size()) thisWire = maybewire[0];

      for (const auto &item : sc->TDCIDEMap()) {
        for (const sim::IDE &ide: item.second) {
          // indexing initializes empty vector
          ret[abs(ide.trackID)].push_back({thisWire, item.first, &ide});
        }
      }
    }
    return ret;
}

std::map<int, std::vector<art::Ptr<recob::Hit>>> sbn::PrepTrueHits(const std::vector<art::Ptr<recob::Hit>> &allHits, 
  const detinfo::DetectorClocksData &clockData, const cheat::BackTrackerService &backtracker) {
  std::map<int, std::vector<art::Ptr<recob::Hit>>> ret;
  for (const art::Ptr<recob::Hit> h: allHits) {
    for (int ID: backtracker.HitToTrackIds(clockData, *h)) {
      ret[abs(ID)].push_back(h);
    }
  }
  return ret;
}

caf::Wall_t sbn::GetWallCross(const geo::BoxBoundedGeo &volume, const TVector3 p0, const TVector3 p1) {
  TVector3 direction = (p1 - p0) * ( 1. / (p1 - p0).Mag());
  std::vector<TVector3> intersections = volume.GetIntersections(p0, direction);

  assert(intersections.size() == 2);

  // get the intersection point closer to p0
  int intersection_i = ((intersections[0] - p0).Mag() < (intersections[1] - p0).Mag()) ? 0 : 1;

  double eps = 1e-3;
  if (abs(intersections[intersection_i].X() - volume.MinX()) < eps) {
    //std::cout << "Left\n";
    return caf::kWallLeft;
  }
  else if (abs(intersections[intersection_i].X() - volume.MaxX()) < eps) {
    //std::cout << "Right\n";
    return caf::kWallRight;
  }
  else if (abs(intersections[intersection_i].Y() - volume.MinY()) < eps) {
    //std::cout << "Bottom\n";
    return caf::kWallBottom;
  }
  else if (abs(intersections[intersection_i].Y() - volume.MaxY()) < eps) {
    //std::cout << "Top\n";
    return caf::kWallTop;
  }
  else if (abs(intersections[intersection_i].Z() - volume.MinZ()) < eps) {
    //std::cout << "Front\n";
    return caf::kWallFront;
  }
  else if (abs(intersections[intersection_i].Z() - volume.MaxZ()) < eps) {
    //std::cout << "Back\n";
    return caf::kWallBack;
  }
  else assert(false);
  //std::cout << "None\n";

  return caf::kWallNone;
}//GetWallCross

//------------------------------------------
caf::g4_process_ sbn::GetG4ProcessID(const std::string &process_name) {
#define MATCH_PROCESS(name) if (process_name == #name) {return caf::kG4 ## name;}
#define MATCH_PROCESS_NAMED(strname, id) if (process_name == #strname) {return caf::kG4 ## id;}
  MATCH_PROCESS(primary)
  MATCH_PROCESS(CoupledTransportation)
  MATCH_PROCESS(FastScintillation)
  MATCH_PROCESS(Decay)
  MATCH_PROCESS(anti_neutronInelastic)
  MATCH_PROCESS(neutronInelastic)
  MATCH_PROCESS(anti_protonInelastic)
  MATCH_PROCESS(protonInelastic)
  MATCH_PROCESS(hadInelastic)
  MATCH_PROCESS_NAMED(kaon+Inelastic, kaonpInelastic)
  MATCH_PROCESS_NAMED(kaon-Inelastic, kaonmInelastic)
  MATCH_PROCESS_NAMED(kaon+Inelastic, kaonpInelastic)
  MATCH_PROCESS_NAMED(kaon-Inelastic, kaonmInelastic)
  MATCH_PROCESS_NAMED(sigma+Inelastic, sigmapInelastic)
  MATCH_PROCESS_NAMED(sigma-Inelastic, sigmamInelastic)
  MATCH_PROCESS_NAMED(pi+Inelastic, pipInelastic)
  MATCH_PROCESS_NAMED(pi-Inelastic, pimInelastic)
  MATCH_PROCESS_NAMED(xi+Inelastic, xipInelastic)
  MATCH_PROCESS_NAMED(xi-Inelastic, ximInelastic)
  MATCH_PROCESS(kaon0LInelastic)
  MATCH_PROCESS(kaon0SInelastic)
  MATCH_PROCESS(lambdaInelastic)
  MATCH_PROCESS_NAMED(anti-lambdaInelastic, anti_lambdaInelastic)
  MATCH_PROCESS(He3Inelastic)
  MATCH_PROCESS(ionInelastic)
  MATCH_PROCESS(xi0Inelastic)
  MATCH_PROCESS(alphaInelastic)
  MATCH_PROCESS(tInelastic)
  MATCH_PROCESS(dInelastic)
  MATCH_PROCESS(anti_neutronElastic)
  MATCH_PROCESS(neutronElastic)
  MATCH_PROCESS(anti_protonElastic)
  MATCH_PROCESS(protonElastic)
  MATCH_PROCESS(hadElastic)
  MATCH_PROCESS_NAMED(kaon+Elastic, kaonpElastic)
  MATCH_PROCESS_NAMED(kaon-Elastic, kaonmElastic)
  MATCH_PROCESS_NAMED(pi+Elastic, pipElastic)
  MATCH_PROCESS_NAMED(pi-Elastic, pimElastic)
  MATCH_PROCESS(conv)
  MATCH_PROCESS(phot)
  MATCH_PROCESS(annihil)
  MATCH_PROCESS(nCapture)
  MATCH_PROCESS(nKiller)
  MATCH_PROCESS(muMinusCaptureAtRest)
  MATCH_PROCESS(muIoni)
  MATCH_PROCESS(eBrem)
  MATCH_PROCESS(CoulombScat)
  MATCH_PROCESS(hBertiniCaptureAtRest)
  MATCH_PROCESS(hFritiofCaptureAtRest)
  MATCH_PROCESS(photonNuclear)
  MATCH_PROCESS(muonNuclear)
  MATCH_PROCESS(electronNuclear)
  MATCH_PROCESS(positronNuclear)
  MATCH_PROCESS(compt)
  MATCH_PROCESS(eIoni)
  MATCH_PROCESS(muBrems)
  MATCH_PROCESS(hIoni)
  MATCH_PROCESS(ionIoni)
  MATCH_PROCESS(hBrems)
  MATCH_PROCESS(muPairProd)
  MATCH_PROCESS(hPairProd)
  MATCH_PROCESS(LArVoxelReadoutScoringProcess)
  MATCH_PROCESS(Transportation)
  MATCH_PROCESS(msc)
  MATCH_PROCESS(StepLimiter)
  MATCH_PROCESS(RadioactiveDecayBase)
  std::cerr << "Error: Process name with no match (" << process_name << ")\n";
  assert(false);
  return caf::kG4UNKNOWN; // unreachable in debug mode
#undef MATCH_PROCESS
#undef MATCH_PROCESS_NAMED
}//GetG4ProcessID
