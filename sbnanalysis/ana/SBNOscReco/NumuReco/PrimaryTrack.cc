#include "PrimaryTrack.h"

int numu::SelectMuon(const std::map<size_t, numu::RecoTrack> &tracks, const numu::RecoInteraction &interaction) {
  if (interaction.slice.primary_index < 0) return -1;

  int muon_candidate = -1;
  float muon_candidate_length = -1;
  
  const numu::RecoParticle &neutrino = interaction.slice.particles.at(interaction.slice.primary_index);
  for (unsigned ID: interaction.PrimaryTracks(tracks)) {
    const numu::RecoTrack &track = tracks.at(ID);
    if (track.pdgid == 13 && track.length > muon_candidate_length) {
      muon_candidate = track.ID;
    }
  }
  return muon_candidate;
}

int numu::SelectLongestTrack(const std::map<size_t, numu::RecoTrack> &tracks, const numu::RecoInteraction &interaction) {
  if (interaction.slice.primary_index < 0) return -1;

  const numu::RecoParticle &neutrino = interaction.slice.particles.at(interaction.slice.primary_index);
  double max_len = -1.;
  int ret = -1;
  for (unsigned ID: interaction.PrimaryTracks(tracks)) {
    const numu::RecoTrack &track = tracks.at(ID);
    if (ret < 0 || track.length > max_len) {
      max_len = track.length;
      ret = track.ID;
    }
  }

  return ret;
}

int numu::SelectLongestIDdMuon(const std::map<size_t, numu::RecoTrack> &tracks, const numu::RecoInteraction &interaction) {
  if (interaction.slice.primary_index < 0) return -1;

  const numu::RecoParticle &neutrino = interaction.slice.particles.at(interaction.slice.primary_index);
  double max_len = -1.;
  int ret = -1;
  for (unsigned ID: interaction.PrimaryTracks(tracks)) {
    const numu::RecoTrack &track = tracks.at(ID);
    // muon ID: exiting or chi2 hypothesis smaller than proton
    if (!track.is_contained || track.chi2_muon < track.chi2_proton) {
      if (ret < 0 || track.length > max_len) {
        max_len = track.length;
        ret = track.ID;
      }
    }
  }

  return ret;
}

