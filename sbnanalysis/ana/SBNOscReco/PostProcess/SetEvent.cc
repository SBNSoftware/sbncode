#include "SetEvent.h"
#include "../NumuReco/PrimaryTrack.h"
#include "../NumuReco/TrackAlgo.h"
#include "../NumuReco/TruthMatch.h"

void ana::SBNOsc::SetEvent(numu::RecoEvent &event, const event::Event &core, const ana::SBNOsc::Cuts &cuts, numu::MCType file_type, bool true_particle_id) {
  // update each reco Interaction to a smarter primary track selector
  unsigned i = 0;
  while (i < event.reco.size()) {
    for (size_t particle_ind: event.reco[i].slice.tracks) {
      numu::RecoTrack &track = event.tracks.at(particle_ind);

      // set containment
      if (cuts.InCalorimetricContainment(track.start) && cuts.InCalorimetricContainment(track.end)) {
        track.is_contained = true;
      }
      else {
        track.is_contained = false;
      }
    }

    // set particle ID
    if (true_particle_id) {
      numu::ApplyTrueParticleID(event.reco[i], event.tracks, event.particles);
    }
    else {
      numu::ApplyParticleID(event.reco[i], event.tracks);
    }

    // set primary track
    int primary_track = numu::SelectMuon(event.tracks, event.reco[i]);

    // remove vertices without a good primary track
    if (primary_track < 0) {
      event.reco.erase(event.reco.begin() + i); 
      continue;
    }

    event.reco[i].primary_track_index = primary_track;

    // set number of true particle types
    event.reco[i].npion = 0;
    event.reco[i].nkaon = 0;
    event.reco[i].nproton = 0;
    for (unsigned ID: event.reco[i].PrimaryTracks(event.tracks)) {
      const numu::RecoTrack &track = event.tracks.at(ID);
      // apply the true particle ID

      if (track.pdgid == 211) {
        event.reco[i].npion ++;
      }
      else if (track.pdgid == 321) {
        event.reco[i].nkaon ++;
      }
      else if (track.pdgid == 2212) {
        event.reco[i].nproton ++;
      }
    }

    // if this is an in-time cosmic file, update the cosmic mode
    // if (file_type == numu::fIntimeCosmic) {
    //   assert(event.reco[i].slice.match.mode == numu::mCosmic || event.reco[i].slice.match.mode == numu::mOther);
    //   event.reco[i].slice.match.mode = numu::mIntimeCosmic;
    // }

    i += 1;
  }
}

