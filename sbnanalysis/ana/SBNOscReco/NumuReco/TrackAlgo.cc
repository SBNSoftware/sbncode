#include "TrackAlgo.h"
#include <numeric>


void numu::ApplyTrueParticleID(const numu::RecoInteraction &interaction, std::map<size_t, numu::RecoTrack> &tracks, const std::map<size_t, numu::TrueParticle> &particles) {
  for (unsigned ID: interaction.PrimaryTracks(tracks)) {
    // apply the true particle ID
    numu::RecoTrack &track = tracks.at(ID);

    int pdg = -1; // "NULL" value
    int match_id = track.truth.GetPrimaryMatchID();
    if (particles.count(match_id)) {
      pdg = particles.at(match_id).pdgid;
    }

    track.pdgid = pdg;
  }
}

void ApplyMinHadronChi2(numu::RecoTrack &track) {
  std::vector<std::pair<int, float>> chi2s {{2212, track.chi2_proton}, {321, track.chi2_kaon}, {211, track.chi2_pion}};

  track.pdgid = 
    std::min_element(chi2s.begin(), chi2s.end(),
      [](auto const &a, auto const &b) {
        return a.second < b.second;
      })->first;
}

void numu::ApplyParticleID(const numu::RecoInteraction &interaction, std::map<size_t, numu::RecoTrack> &tracks) {
  // example algorithm borrowed from Rhiannon Jones
  // first look for a track with length > 100cm that exits 
  // call this the muon
  int exiting_muon_candidate = -1;
  for (unsigned ID: interaction.PrimaryTracks(tracks)) {
    const numu::RecoTrack &track = tracks.at(ID);
    if (track.length > 100 && !track.is_contained) {
      if (exiting_muon_candidate != -1) { // two candidates -- don't apply this ID
        exiting_muon_candidate = -1;
        break;
      }
      exiting_muon_candidate = track.ID;
    }
  }

  if (exiting_muon_candidate != -1) {
    for (unsigned ID: interaction.PrimaryTracks(tracks)) {
      numu::RecoTrack &track = tracks.at(ID);
      if (ID == exiting_muon_candidate) {
        track.pdgid = 13; // muon
      }
      else {
        ApplyMinHadronChi2(track);
      }
    }
    // Done!
    return;
  }

  // otherwise:

  std::vector<std::pair<size_t, float>> muon_candidates;
  for (unsigned ID: interaction.PrimaryTracks(tracks)) {
    const numu::RecoTrack &track = tracks.at(ID);
    bool very_long = true;
    // get the subleading max length
    for (unsigned otherID: interaction.PrimaryTracks(tracks)) {
    const numu::RecoTrack &other = tracks.at(otherID);
      if (other.ID == track.ID) continue; // same track
      if (track.length < other.length * 1.5) {
        very_long = false;
        break;
      }
    }
    if (very_long || (track.chi2_muon < 16 && track.chi2_proton > 80)) {
      std::pair<size_t, float> pair {ID, track.chi2_muon};
      muon_candidates.push_back(pair);
    }
  }

  // get the best muon candidate
  int muon = -1;
  if (muon_candidates.size() > 0) {
    muon = std::min_element(muon_candidates.begin(), muon_candidates.end(), 
      [](const auto &a, const auto &b) {
        return a.first < b.first;
      })->first;
  }

  for (unsigned ID: interaction.PrimaryTracks(tracks)) {
    numu::RecoTrack &track = tracks.at(ID);
    if (ID == muon) {
      track.pdgid = 13;
    }
    else if (track.chi2_proton < 80.) {
      track.pdgid = 2212;
    }
    else {
      ApplyMinHadronChi2(track);
    }
  }
  return;
}

// TODO: make this more sophisticated
float numu::TrackMomentum(const numu::RecoTrack &track) {
  if (track.is_contained) {
    return numu::RangeMomentum(track);
  }
  else {
    return numu::MCSMomentum(track);
  }
}

// TODO: make this more sophisticated
float numu::RangeMomentum(const numu::RecoTrack &track) {
  return track.range_momentum_muon;
}

// TODO: make this more sophisticated
float numu::MCSMomentum(const numu::RecoTrack &track) {
  return track.mcs_muon.fwd_mcs_momentum;
}


float numu::MeanTruncateddQdx(const anab::Calorimetry &calo) {
  // copy the dQdx list for partial sorting to find median
  std::vector<float> dQdx = calo.dQdx();

  // if no or only 1 calo points, then no dQdx
  if (dQdx.size() <= 1) return -9999.;

  // calculate the median
  float median = -1;
  if (dQdx.size() % 2 == 0 /* even */) {
    const auto median_lo = dQdx.begin() + dQdx.size()/2 - 1;
    const auto median_hi = dQdx.begin() + dQdx.size()/2;
    std::nth_element(dQdx.begin(), median_lo, dQdx.end());
    float median_lo_val = *median_lo;
    std::nth_element(dQdx.begin(), median_hi, dQdx.end());
    float median_hi_val = *median_hi;
    median = (median_lo_val + median_hi_val) / 2.;    
  }
  else /* odd */ {
    const auto median_ptr = dQdx.begin() + dQdx.size()/2;
    std::nth_element(dQdx.begin(), median_ptr, dQdx.end());
    median = *median_ptr;
  }

  // get the mean and variance
  float mean = std::accumulate(dQdx.begin(), dQdx.end(), 0.) / dQdx.size();  
  float var = 0.;
  for (float d: dQdx) {
    var += (d - mean) * (d - mean);
  }
  var = var / dQdx.size();

  // truncated mean
  float trunc_sum = 0.;
  unsigned n_point = 0;
  for (float d: dQdx) {
    if (var < (d - mean) * (d - mean)) {
      trunc_sum += d;
      n_point ++;
    }
  }
  // if no points available, return garbage
  if (n_point == 0) return -9999.;

  return trunc_sum / n_point;
}
