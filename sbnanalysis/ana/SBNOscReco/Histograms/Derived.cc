#include "Derived.h"

float numu::dist2Match(const event::Interaction &truth, const std::vector<numu::RecoInteraction> &candidates) {
  float dist = -1;
  for (const numu::RecoInteraction &candidate: candidates) {
    float this_dist = (candidate.position - truth.neutrino.position).Mag();
    if (dist < 0 || this_dist < dist) dist = this_dist;
  }
  return dist;

}

float numu::trackMatchCompletion(unsigned g4id, const numu::RecoEvent &event) {
  if (g4id == -1) return -1;

  float completion = -1;
  float most_matched_energy = 0.; 
  float this_energy = event.particles.at(g4id).deposited_energy;
  for (const auto &pair: event.tracks) {
    const numu::RecoTrack &track = pair.second;
    if (track.truth.GetPrimaryMatchID() == g4id) {
      float this_matched_energy = track.truth.matches[0].energy;
      if (this_matched_energy > most_matched_energy) {
        completion = this_matched_energy / this_energy;
        most_matched_energy = this_matched_energy;
      }
    }
  }
  return completion;
}
