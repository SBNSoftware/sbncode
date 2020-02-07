#include "TruthMatch.h"
#include "nusimdata/SimulationBase/MCTruth.h"

struct MatchSet {
  std::vector<bool> match_is_primary;
  std::vector<int> reco_to_truth;
};

static const bool _verbose_truth_match = true;

// declare heuristics
MatchSet ApplyNumuCCHeuristic(const numu::RecoEvent &event, const event::Interaction &truth, unsigned truth_i); 

// Method of truth matching:
//
// First iterate through the true neutrino interactions and use some heuristics to
// see if there are any obvious matches between the true neutrino interaction and
// the reco slices. 
//
// The current heuristics in use are:
//   1. numu CC: Match the reco slice that best reconstructions the final state muon 
//   2. TODO: other heuristics (nue CC? NC?)
//
// After applying these heuristics, iterate through the reco interactions and match
// each reconstructed interaction to the truth interaction that it matches the most.
//
// We don't allow multiple reconstructed interactions to match the same true
// neutrino interaction. In the case that multiple reco interactions match to the
// same true neutrino, we break the tie with the following algorithm:
//
// Tiebreakers:
//   1. If one of the reco interactions satisfied one of the heuristics, then that
//   interaction matches the true neutrino (the primacy of each heurestic is in the
//   order enumerated above)
//   2. If none of the reco interactions satisfy a heuristic, then we match the one
//   that matches the most __primary__ energy of the neutrino interaction. We only
//   consider primary energy to remove reco interactions that are associated with
//   downstream energy of the main reco neutrino interaction, such as a neutron from
//   the neutrino vertex that progogates in the active volume and re-interacts
//   inelastically. 
//
// In the case of a tiebreaker being applied, all reco interactions get a
// "NonPrimary" interaction mode (i.e. CCNonPrimary, NCNonPrimary). A reco
// interaction will also be assigned a "NonPrimary" interaction mode if it is
// primarily comprised of energy from non-primary particles (see the motivating case
// above)
MatchSet FindMatches(const numu::RecoEvent &event, const std::vector<event::Interaction> &truth) {
  // matching from reco -> truth (-1 means cosmic)
  std::vector<int> reco_to_truth(event.reco.size(), -1);
  std::vector<bool> match_is_primary(event.reco.size(), false);

  // matching using heuristics
  std::vector<int> reco_to_truth_h(event.reco.size(), -1);
  // if heuristic is primary
  std::vector<bool> match_is_primary_h(event.reco.size(), false);

  // matching using energy
  std::vector<int> truth_best_energy_matched(truth.size(), 0.);

  if (_verbose_truth_match) std::cout << "\n\nNew truth matching!\n";
  if (_verbose_truth_match) std::cout << "NRECO: " << event.reco.size() << " NTRUTH: " << truth.size() << std::endl; 

  // find the best match for each true interaction
  for (unsigned truth_i = 0; truth_i < truth.size(); truth_i++) {
    const event::Interaction &t = truth[truth_i];

    if (_verbose_truth_match) std::cout << "Trying heuristic for neutrino interaction: " << truth_i << std::endl;

    MatchSet matches; 

    // Hueristic for numu CC
    if (t.neutrino.iscc && abs(t.neutrino.pdg) == 14) {
      matches = ApplyNumuCCHeuristic(event, t, truth_i);
    }
    // TODO: other heuristics???
    else {}

    // merge heuristics
    if (matches.reco_to_truth.size() > 0)  {
      for (unsigned reco_i = 0; reco_i < event.reco.size(); reco_i++) {
        if (matches.reco_to_truth[reco_i] != -1) {
          assert(reco_to_truth_h[reco_i] == -1);
          reco_to_truth_h[reco_i] = matches.reco_to_truth[reco_i];
          match_is_primary_h[reco_i] = matches.match_is_primary[reco_i];
        }
      }
    }
  }

  // calculate the true primary deposited energy for each true interaction
  std::vector<float> true_primary_energy(truth.size(), 0.);
  for (unsigned truth_i = 0; truth_i < truth.size(); truth_i++) {
    for (const event::FinalStateParticle &p: truth[truth_i].finalstate) {
      if (p.is_primary && event.particles.count(p.G4ID)) {
        true_primary_energy[truth_i] += event.particles.at(p.G4ID).deposited_energy;
      }
    }
  }

  // confirm heuristic matches
  for (unsigned reco_i = 0; reco_i < event.reco.size(); reco_i++) {
    if (reco_to_truth_h[reco_i] != -1) {
      reco_to_truth[reco_i] = reco_to_truth_h[reco_i];
      match_is_primary[reco_i] = match_is_primary_h[reco_i];
    }
  }

  // Post heuristics: for the interactions without a match already,
  // match using deposited energy
  for (unsigned reco_i = 0; reco_i < event.reco.size(); reco_i++) {
    const numu::RecoInteraction &reco = event.reco[reco_i];
    // match already made
    if (reco_to_truth_h[reco_i] != -1)  continue;

    // find the total energy matching each true interaction
    std::vector<float> primary_matching_energy(truth.size(), 0.);
    std::vector<float> matching_energy(truth.size(), 0.);

    float this_total_deposited_energy = 0.;
    // TODO: this only includes tracks right now -- also include showers?
    for (size_t trackID: event.reco[reco_i].slice.tracks) {
      const numu::RecoTrack &track = event.tracks.at(trackID);
      int particle_match_id = track.truth.GetPrimaryMatchID();
      if (particle_match_id >= 0) {
        const numu::TrueParticle &particle = event.particles.at(particle_match_id);
        if (particle.interaction_id >= 0) {
          float matched_deposited_energy = track.truth.matches[0].energy;
          matching_energy[particle.interaction_id] += matched_deposited_energy;
          if (particle.IsPrimary()) primary_matching_energy[particle.interaction_id] += matched_deposited_energy;
        } 
      } 
      this_total_deposited_energy += track.truth.total_deposited_energy; 
    }

    // find if any match accounts for __most__ of the energy
    // if none does, then this interaction is matched to a cosmic
    int candidate_match = -1;
    float neutrino_matched_energy = 0.;
    for (unsigned truth_i = 0; truth_i < truth.size(); truth_i++) {
       if (_verbose_truth_match) std::cout << "Reco: " << reco_i << " to truth: " << truth_i << " frac match: " << (matching_energy[truth_i] / this_total_deposited_energy) << std::endl;

      if (matching_energy[truth_i] / this_total_deposited_energy > 0.5) {
         if (_verbose_truth_match) std::cout << "Reco interaction: " << reco_i << " match truth interaction: " << truth_i << " with match fraction: " << (matching_energy[truth_i] / this_total_deposited_energy) << " primary frac: " << (primary_matching_energy[truth_i] / matching_energy[truth_i]) << std::endl; 

        // most energy matched to non-primary: means this match is non-primary
        if (primary_matching_energy[truth_i] / matching_energy[truth_i] < 0.5) {
          if (_verbose_truth_match) std::cout << "non-primary match\n";
          reco_to_truth[reco_i] = truth_i; 
          match_is_primary[reco_i] = false;
          continue;
        }

        // Attempting primary match -- see if there already was one
        bool is_primary = true; // primary until proven otherwise
        for (unsigned other_reco_i = 0; other_reco_i < event.reco.size(); other_reco_i++) {
          if (reco_to_truth[other_reco_i] == truth_i && match_is_primary[other_reco_i]) {
            if (_verbose_truth_match) std::cout << "Duplicate match! This: " << reco_i << " other: " << other_reco_i << std::endl;

            // if the other one was matched via a heuristic, this one loses
            if (reco_to_truth_h[other_reco_i] != -1) {
              if (_verbose_truth_match) std::cout << "Other is heuristic-matched.\n";
              is_primary = false;
              break;
            } 
            // otherwise, see who has more matching primary energy
            //
            // The previous match does
            else if (primary_matching_energy[truth_i] < truth_best_energy_matched[truth_i]) {
              if (_verbose_truth_match) std::cout << "Other has more matching primary energy. This: " << primary_matching_energy[truth_i] << " other: " << truth_best_energy_matched[truth_i] << std::endl;
              is_primary = false;
              break;
            } 
            // this one does! Set the other to secondary
            else {
              if (_verbose_truth_match) std::cout << "This has more matching primary energy. This: " << primary_matching_energy[truth_i] << " other: " << truth_best_energy_matched[truth_i] << std::endl;
              match_is_primary[other_reco_i] = false;
            }
          }
        }
        if (is_primary) {
          truth_best_energy_matched[truth_i] = primary_matching_energy[truth_i];
        }
        // Found the match!
        reco_to_truth[reco_i] = truth_i;
        match_is_primary[reco_i] = is_primary;
        break;
      }
    } 
  } 
  MatchSet s;
  s.reco_to_truth = std::move(reco_to_truth);
  s.match_is_primary = std::move(match_is_primary);
  return std::move(s);
}

numu::SliceTruth MakeSliceTruth(const numu::RecoEvent &event, const std::vector<event::Interaction> &truth, unsigned reco_i, int interaction_id, bool match_is_primary) {
  numu::SliceTruth ret;
  ret.interaction_id = interaction_id;

  if (_verbose_truth_match) std::cout << "Reco: " << reco_i << " matched to: ";
  if (interaction_id != -1) {
    assert(event.type == numu::fOverlay);
    const event::Interaction &t = truth[interaction_id];
    if (t.neutrino.iscc && match_is_primary) {
      if (_verbose_truth_match) std::cout << "CC";
      ret.mode = numu::mCC;
    }
    else if (t.neutrino.iscc && !match_is_primary) {
      if (_verbose_truth_match) std::cout << "CC (non primary)";
      ret.mode = numu::mCCNonPrimary;
    }
    else if (!t.neutrino.iscc && match_is_primary) {
      if (_verbose_truth_match) std::cout << "NC";
      ret.mode = numu::mNC;
    }
    else if (!t.neutrino.iscc && !match_is_primary) {
      if (_verbose_truth_match) std::cout << "NC (non primary)";
      ret.mode = numu::mNCNonPrimary;
    }
    else assert(false); // unreachable
  }
  else if (event.type == numu::fIntimeCosmic) {
    if (_verbose_truth_match) std::cout << "Intime cosmic"; 
    ret.mode = numu::mIntimeCosmic; 
  }  
  else {
    if (_verbose_truth_match) std::cout << "Cosmic"; 
    ret.mode = numu::mCosmic;
  }
  std::cout << std::endl;

  float total_deposited_energy = 0.;
  float matching_neutrino_energy = 0.;
  float cosmic_energy = 0.;
  float unmatched_energy = 0.;
  // calculate matching energy
  for (size_t ind: event.reco[reco_i].slice.tracks) {
    const numu::RecoTrack &track = event.tracks.at(ind); 
    float this_matched_energy = 0.;
    for (const numu::TrackTruth::ParticleMatch &pmatch: track.truth.matches) {
      this_matched_energy += pmatch.energy;

      // FIXME: Sometimes particles that have true energy depositions in the detector are
      //        not saved by G4 (I am not sure why). This seems to not happen for the 
      //        primary match to a given track, so we can usually assert that the 
      //        particle matches for the primary match. However, sometimes particles
      //        with small energy depositions seem to be thrown away. In this case,
      //        we cannot determine the origin of the energy, so we count it 
      //        as "unmatched" 
      if (event.particles.count(pmatch.G4ID)) {
        const numu::TrueParticle &particle = event.particles.at(pmatch.G4ID);
        if (particle.interaction_id < 0.) {
          cosmic_energy += pmatch.energy;
        }
        else if (particle.interaction_id == interaction_id) {
          matching_neutrino_energy += pmatch.energy;
        }
      }
    }
    unmatched_energy = track.truth.total_deposited_energy - this_matched_energy;
    total_deposited_energy += track.truth.total_deposited_energy;
  }

  if (_verbose_truth_match) std::cout << "Cosmic match frac: " << (cosmic_energy/total_deposited_energy) << " Neutrino frac: " << (matching_neutrino_energy/total_deposited_energy) << " Unmatched frac: " << (unmatched_energy/total_deposited_energy) << std::endl;

  ret.total_deposited_energy = total_deposited_energy;
  ret.neutrino_match_energy = matching_neutrino_energy;
  ret.cosmic_energy = cosmic_energy;

  ret.unmatched_energy = unmatched_energy;
  // Because an unmatched energy of zero requires different floating
  // point values to cancel to zero, sometimes floating point error
  // can create a small artifical  "unmatched_energy" value (which
  // can be negative)
  //
  // To avoid confusion, round very small values of unmatched_energy
  // to 0 
  if (abs(ret.unmatched_energy / ret.total_deposited_energy) /* make unit-less */ < 1e-5) {
    ret.unmatched_energy = 0.;
  }

  return ret;
} 

void numu::ApplySliceTruthMatch(numu::RecoEvent &event, const std::vector<event::Interaction> &truth) {
  MatchSet matches = FindMatches(event, truth);

  // coherence check
  std::vector<bool> has_primary_match(truth.size(), false);
  for (unsigned i = 0; i < event.reco.size(); i++) {
    if (matches.reco_to_truth[i] != -1 && matches.match_is_primary[i]) {
      assert(!has_primary_match[matches.reco_to_truth[i]]);    
      has_primary_match[matches.reco_to_truth[i]] = true;
    }
  }

  // now apply the matches and calculate the truth
  for (unsigned reco_i = 0; reco_i < event.reco.size(); reco_i++) {
    event.reco[reco_i].slice.truth = MakeSliceTruth(event, truth, reco_i, matches.reco_to_truth[reco_i], matches.match_is_primary[reco_i]); 
  }
}


// heuristics
MatchSet ApplyNumuCCHeuristic(const numu::RecoEvent &event, const event::Interaction &truth, unsigned truth_i) {
  // matching using heuristics
  std::vector<int> reco_to_truth_h(event.reco.size(), -1);
  // if heuristic is primary
  std::vector<bool> match_is_primary_h(event.reco.size(), false);

  if (_verbose_truth_match) std::cout << "numu CC heuristic matches\n";

  // get the reco candidates -- look for any match to the muon
  std::vector<int> muon_matches(event.reco.size(), -1);
  for (unsigned reco_i = 0; reco_i < event.reco.size(); reco_i++) {
    const numu::RecoInteraction &reco = event.reco[reco_i];
    const numu::RecoParticle &neutrino = reco.slice.particles.at(reco.slice.primary_index);
    std::cout << "Daughters of reco: (" << reco_i << ") ";
    for (size_t ind: neutrino.daughters) {
      std::cout << ind << " ";
      if (event.tracks.count(ind)) {
        float dist = (reco.position - event.tracks.at(ind).start).Mag();
        std::cout << "(" << dist << ") ";
      }
    }
    std::cout << std::endl;
    for (unsigned ID: event.reco[reco_i].PrimaryTracks(event.tracks)) {
      const numu::RecoTrack &track = event.tracks.at(ID);

      if (track.truth.GetPrimaryMatchID() == truth.lepton.G4ID) {
        float energy_match_frac = track.truth.matches[0].energy / event.particles.at(truth.lepton.G4ID).deposited_energy; 
        if (_verbose_truth_match) {
          std::cout << "heuristic check for reco: " << reco_i << ". True muon energy: " << event.particles.at(truth.lepton.G4ID).deposited_energy << " deposited: " << event.particles.at(truth.lepton.G4ID).energy_loss << " match energy: " << track.truth.matches[0].energy << 
	      ". True muon start X: " << event.particles.at(truth.lepton.G4ID).start.X() << 
	      ". True muon start Y: " << event.particles.at(truth.lepton.G4ID).start.Y() << 
	      ". True muon start Z: " << event.particles.at(truth.lepton.G4ID).start.Z() << 
	      " match start X: " << track.start.X() <<
	      " match start Y: " << track.start.Y() <<
	      " match start Z: " << track.start.Z() <<
	      " match end X: " << track.end.X() <<
	      " match end Y: " << track.end.Y() <<
	      " match end Z: " << track.end.Z() << 
              " match length: " << track.length << 
              " match ID: " << track.ID << std::endl;
        }
        if (energy_match_frac > 0.05) {

          // muon got split in this reco -- choose the one closer to the start
          if (muon_matches[reco_i] != -1) {
            std::cout << "BAD: Two reco tracks matching muon id: " << truth.lepton.G4ID << std::endl;
            const numu::RecoTrack &this_track = track;
            std::cout << "this track: " << "length: " << this_track.length << " match energy frac: " << (this_track.truth.matches[0].energy / event.particles.at(truth.lepton.G4ID).deposited_energy) << " x: " << this_track.start.X() << " y: " << this_track.start.Y() << " z: " << this_track.start.Z() << std::endl; 
            const numu::RecoTrack &other_track = event.tracks.at(muon_matches[reco_i]);
            std::cout << "other track: " << "length: " << other_track.length << " match energy frac: " << (other_track.truth.matches[0].energy / event.particles.at(truth.lepton.G4ID).deposited_energy) << " x: " << other_track.start.X() << " y: " << other_track.start.Y() << " z: " << other_track.start.Z() << std::endl; 

            float this_distance = std::min((event.particles.at(truth.lepton.G4ID).start - this_track.start).Mag(),
                                           (event.particles.at(truth.lepton.G4ID).start - this_track.end).Mag());
            float other_distance = std::min((event.particles.at(truth.lepton.G4ID).start - other_track.start).Mag(),
                                            (event.particles.at(truth.lepton.G4ID).start - other_track.end).Mag());
            if (this_distance < other_distance) {
              muon_matches[reco_i] = track.ID;
            }
          }
          // not split (yet)
          else {
            muon_matches[reco_i] = track.ID;
          }
        }
        else if (_verbose_truth_match) std::cout << "HEURISTIC: match energy too low\n";
      }
    } 
  }

  unsigned nmatch = 0;
  for (int match: muon_matches) {
    if (match != -1) nmatch ++;
  }

  if (_verbose_truth_match) std::cout << "matches for heuristic: " << nmatch << std::endl;

  // one match -- great! This one is the reco-truth match
  if (nmatch == 1) {
    for (unsigned reco_i = 0; reco_i < event.reco.size(); reco_i++) {
      if (muon_matches[reco_i] != -1) {
        assert(reco_to_truth_h[reco_i] == -1);
        reco_to_truth_h[reco_i] = truth_i;
        match_is_primary_h[reco_i] = true;
     }
    }
  }
  // multiple matches -- use tiebreaker
  //
  // If there are multiple tracks matching to the muon then
  // it most likely got split in two. Tiebreak by using the track
  // with the start point closer to the true vertex. The loser of
  // the tie-breaker is set as non-primary
  else if (nmatch > 1) {
    float distance_match = -1;
    int best_index = -1;
    for (unsigned reco_i = 0; reco_i < event.reco.size(); reco_i++) {
      int match = muon_matches[reco_i];
      if (match != -1) {
        float this_distance = std::min((event.tracks.at(match).start - truth.neutrino.position).Mag(),
                                       (event.tracks.at(match).end   - truth.neutrino.position).Mag());

        if (_verbose_truth_match) std::cout << "reco match: " << reco_i << " has distance from truth: " << this_distance << std::endl;

        if (distance_match < 0. || this_distance < distance_match) {
          distance_match = this_distance; 
          best_index = reco_i;
        }
      }
    }
    assert(best_index != -1);
    if (_verbose_truth_match) std::cout << "best reco match to heuristic is: " << best_index << std::endl;
    for (unsigned reco_i = 0; reco_i < event.reco.size(); reco_i++) {
      int match = muon_matches[reco_i];
      if (match != -1) {
        reco_to_truth_h[reco_i] = truth_i;
        match_is_primary_h[reco_i] = best_index == reco_i;
      }
    }
  }
  // no-match :(
  else {}

  return {match_is_primary_h, reco_to_truth_h};
}


