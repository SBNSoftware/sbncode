#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Hit.h"

#include "TrackReducer.h"
#include "Data.h"

#include "ana/SBNOscReco/RecoUtils/RecoUtils.h"

void sbnana::TrackReducer::Initialize(fhicl::ParameterSet* config) {
  fTree->Branch("tracks", &fTracks);
}

bool sbnana::TrackReducer::ProcessEvent(const gallery::Event& ev, const std::vector<event::Interaction> &truth, std::vector<event::RecoInteraction>& reco) {
  std::cout << "New Event!\n";
  fTracks.reco_tracks.clear();
  fTracks.reco_showers.clear();
  fTracks.truth.clear();

  // map of particle track ID's to indices
  std::map<int, unsigned> trackID_to_index;

  // get the important truth stuff
  const std::vector<simb::MCParticle> mcparticles = *ev.getValidHandle<std::vector<simb::MCParticle>>("largeant");

  unsigned ind = 0;
  for (const simb::MCParticle &mc_particle: mcparticles) {
    trackID_to_index[mc_particle.TrackId()] = ind;
    ind ++; 
  }

  for (const simb::MCParticle &mc_particle: mcparticles) {
    TrueParticle particle;
    particle.pid = mc_particle.PdgCode();
    particle.is_primary = (mc_particle.Process() == "primary");
    particle.is_michel = abs(mc_particle.PdgCode()) == 11 && (mc_particle.Process() == "Decay" || mc_particle.Process() == "muMinusCaptureAtRest");
    particle.is_delta  = abs(mc_particle.PdgCode()) == 11 && (mc_particle.Process() == "muIoni");
    for (unsigned i = 0; i < mc_particle.NumberDaughters(); i++) {
      int daughter_id = mc_particle.Daughter(i);
      if (trackID_to_index.count(daughter_id)) particle.daughters.push_back(trackID_to_index.at(daughter_id));
    }
    particle.energy = mc_particle.E();

    if (mc_particle.NumberTrajectoryPoints() > 0) {
      for (unsigned i = 0; i < mc_particle.NumberTrajectoryPoints(); i++) {
        particle.trajectory.push_back(mc_particle.Position(i).Vect());
      }
    }
    else particle.trajectory.push_back(mc_particle.Position().Vect());

    std::cout << "Particle: " << mc_particle.TrackId() << " PID: " << mc_particle.PdgCode() << " process: " << mc_particle.Process() << " energy: " << mc_particle.E() << std::endl;
    fTracks.truth.push_back(std::move(particle));

  }

  const auto &track_handle = ev.getValidHandle<std::vector<recob::Track>>("pandoraTrack");
  const auto &shower_handle = ev.getValidHandle<std::vector<recob::Shower>>("pandoraShower");

  art::FindManyP<recob::Hit> tracks_to_hits(track_handle, ev, "pandoraTrack");
  art::FindManyP<recob::PFParticle> tracks_to_particles(track_handle, ev, "pandoraTrack");

  art::FindManyP<recob::PFParticle> showers_to_particles(shower_handle, ev, "pandoraShower");
  art::FindManyP<recob::Hit> showers_to_hits(shower_handle, ev, "pandoraShower");

  // setup ID mappings
  std::map<unsigned, unsigned> particle_to_showerID;
  std::map<unsigned, unsigned> particle_to_trackID;
  for (unsigned track_id = 0; track_id < track_handle->size(); track_id++) {
    unsigned particleID = tracks_to_particles.at(track_id).at(0)->Self();
    std::cout << "track: " << track_id << " to particle: " << particleID << std::endl;
    particle_to_trackID[particleID] = track_id;
  }

  for (unsigned shower_id = 0; shower_id < shower_handle->size(); shower_id++) {
    unsigned particleID = showers_to_particles.at(shower_id).at(0)->Self();
    std::cout << "shower: " << shower_id << " to particle: " << particleID << std::endl;
    particle_to_showerID[particleID] = shower_id;
  }

  for (unsigned track_id = 0; track_id < track_handle->size(); track_id++) {
    const recob::Track &rb_track = (*track_handle)[track_id];
    Track track;
    for (unsigned i = 0; i < rb_track.CountValidPoints(); i++) {
      TVector3 traj (rb_track.LocationAtPoint(i).X(), rb_track.LocationAtPoint(i).Y(), rb_track.LocationAtPoint(i).Z());
      track.trajectory.push_back(traj);
    }
    std::vector<art::Ptr<recob::Hit>> hits = tracks_to_hits.at(track_id);

    // truth match
    int mcp_track_id =  SBNRecoUtils::TrueParticleIDFromTotalTrueEnergy(*fProviderManager, hits, true); 
    if (trackID_to_index.count(mcp_track_id)) track.truth_match = trackID_to_index.at(mcp_track_id);
    else track.truth_match = -1;

    // get daughters
    const recob::PFParticle &pfparticle = *tracks_to_particles.at(track_id).at(0); 
    for (unsigned daughter: pfparticle.Daughters()) {
      std::cout << "Daughter! " << daughter << std::endl;
      if (particle_to_trackID.count(daughter)) track.track_daughters.push_back(particle_to_trackID.at(daughter)); 
      else if (particle_to_showerID.count(daughter)) track.shower_daughters.push_back(particle_to_showerID.at(daughter)); 
    }

    std::cout << "New Track! n traj: " << track.trajectory.size() << std::endl;;

    fTracks.reco_tracks.push_back(track);
  }

  for (unsigned shower_id = 0; shower_id < shower_handle->size(); shower_id++) {
    const recob::Shower &rb_shower = (*shower_handle)[shower_id];
    Shower shower;
    shower.direction = rb_shower.Direction();
    shower.direction_err = rb_shower.DirectionErr();
    shower.start = rb_shower.ShowerStart();
    shower.start_err = rb_shower.ShowerStartErr();

    std::vector<art::Ptr<recob::Hit>> hits = showers_to_hits.at(shower_id);

    // truth match
    int mcp_track_id =  SBNRecoUtils::TrueParticleIDFromTotalTrueEnergy(*fProviderManager, hits, true); 
    if (trackID_to_index.count(mcp_track_id)) shower.truth_match = trackID_to_index.at(mcp_track_id);
    else shower.truth_match = -1;

    // get daughters
    const recob::PFParticle &pfparticle = *showers_to_particles.at(shower_id).at(0); 
    for (unsigned daughter: pfparticle.Daughters()) {
      std::cout << "Daughter! " << daughter << std::endl;
      if (particle_to_trackID.count(daughter)) shower.track_daughters.push_back(particle_to_trackID.at(daughter)); 
      else if (particle_to_showerID.count(daughter)) shower.shower_daughters.push_back(particle_to_showerID.at(daughter)); 
    }

    fTracks.reco_showers.push_back(shower);
    std::cout << "New Shower! at: " << rb_shower.ShowerStart().X() << std::endl;
    std::cout << "Match to: " << mcp_track_id << std::endl;
  }

  return true;
}

DECLARE_SBN_PROCESSOR(sbnana::TrackReducer)
