#include "RecoEvent.h"

std::vector<size_t> numu::RecoInteraction::PrimaryTracks(const std::map<size_t, numu::RecoTrack> &tracks) const {
  std::vector<size_t> ret;

  const numu::RecoParticle &neutrino = slice.particles.at(slice.primary_index);

  for (size_t id: neutrino.daughters) {
    if (!tracks.count(id)) continue;
    float dist = (position - tracks.at(id).start).Mag();
    if (dist > 10.) continue;
    ret.push_back(id);
  }

  return ret;
}

