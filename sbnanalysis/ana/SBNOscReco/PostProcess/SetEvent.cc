#include "SetEvent.h"
#include "../NumuReco/PrimaryTrack.h"
#include "../NumuReco/TruthMatch.h"

void ana::SBNOsc::SetEvent(numu::RecoEvent &event, const event::Event &core, const ana::SBNOsc::Cuts &cuts, numu::MCType file_type, bool use_calorimetry) {
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

    int primary_track;
    if (use_calorimetry) {
      primary_track = numu::SelectLongestIDdMuon(event.tracks, event.reco[i].slice);
    }
    else {
      primary_track = numu::SelectLongestTrack(event.tracks, event.reco[i].slice);
    }

    // remove vertices without a good primary track
    if (primary_track < 0) {
      event.reco.erase(event.reco.begin() + i); 
      continue;
    }

    event.reco[i].primary_track_index = primary_track;

    // if this is an in-time cosmic file, update the cosmic mode
    if (file_type == numu::fIntimeCosmic) {
      assert(event.reco[i].slice.match.mode == numu::mCosmic || event.reco[i].slice.match.mode == numu::mOther);
      event.reco[i].slice.match.mode = numu::mIntimeCosmic;
    }

    i += 1;
  }

  // re-do truth matching
  numu::ApplyPrimaryTrackTruthMatch(event, core.truth);
}

