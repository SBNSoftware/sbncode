#include "Derived.h"

float numu::dist2Match(const event::Interaction &truth, const std::vector<numu::RecoInteraction> &candidates) {
  float dist = -1;
  for (const numu::RecoInteraction &candidate: candidates) {
    float this_dist = (candidate.position - truth.neutrino.position).Mag();
    if (dist < 0 || this_dist < dist) dist = this_dist;
  }
  return dist;

}

float numu::trackMatchCompletion(unsigned truth_index, const numu::RecoEvent &event) {
  for (const numu::RecoInteraction &reco: event.reco) {
    const numu::RecoTrack &primary_track = event.tracks.at(reco.primary_track_index);
    if (reco.slice.match.mctruth_track_id == truth_index && primary_track.match.is_primary) {
      return primary_track.match.completion; 
    }
  }
  return -1;
}
