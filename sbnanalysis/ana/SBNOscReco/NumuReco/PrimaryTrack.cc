#include "PrimaryTrack.h"

int numu::SelectMuon(const std::map<size_t, numu::RecoTrack> &tracks, const numu::RecoSlice &slice) {
  if (slice.primary_index < 0) return -1;

  const numu::RecoParticle &neutrino = slice.particles.at(slice.primary_index);
  for (size_t pfp_index: neutrino.daughters) {
    if (!tracks.count(pfp_index)) continue;
    if (tracks.at(pfp_index).pdgid == 13) return pfp_index;  
  }
  return -1;
}

int numu::SelectLongestTrack(const std::map<size_t, numu::RecoTrack> &tracks, const numu::RecoSlice &slice) {
  if (slice.primary_index < 0) return -1;

  const numu::RecoParticle &neutrino = slice.particles.at(slice.primary_index);
  double max_len = -1.;
  int ret = -1;
  for (size_t pfp_index: neutrino.daughters) {
    const numu::RecoParticle &daughter = slice.particles.at(pfp_index);
    if (tracks.count(daughter.ID)) {
      if (ret < 0 || tracks.at(daughter.ID).length > max_len) {
        max_len = tracks.at(daughter.ID).length;
        ret = daughter.ID;
      }
    }
  }

  return ret;
}

int numu::SelectLongestIDdMuon(const std::map<size_t, numu::RecoTrack> &tracks, const numu::RecoSlice &slice) {
  if (slice.primary_index < 0) return -1;

  const numu::RecoParticle &neutrino = slice.particles.at(slice.primary_index);
  double max_len = -1.;
  int ret = -1;
  for (size_t pfp_index: neutrino.daughters) {
    const numu::RecoParticle &daughter = slice.particles.at(pfp_index);
    if (tracks.count(daughter.ID)) {
      const numu::RecoTrack &track = tracks.at(daughter.ID);
      // muon ID: exiting or chi2 hypothesis smaller than proton
      if (!track.is_contained || track.chi2_muon < track.chi2_proton) {
        if (ret < 0 || track.length > max_len) {
          max_len = track.length;
          ret = daughter.ID;
        }
      }
    }
  }

  return ret;
}

