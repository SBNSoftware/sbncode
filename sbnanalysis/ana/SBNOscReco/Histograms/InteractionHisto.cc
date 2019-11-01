#include "InteractionHisto.h"
#include "Derived.h"

#include "TH1D.h"
#include "TH2D.h"

namespace ana {
  namespace SBNOsc {

void InteractionHistos::Initialize(const std::string &postfix) {
#define INT_HISTO(name, n_bins, lo, hi)    name = new TH1D((#name"_" + postfix).c_str(), #name, n_bins, lo, hi); fAllHistos.push_back(name)


  INT_HISTO(track_length, 101, -10, 1000);
  INT_HISTO(track_p, 50, 0., 5.);
  INT_HISTO(true_deposited_energy, 50., 0., 5.);
  INT_HISTO(nuE, 50, 0., 5.);
  INT_HISTO(beam_center_distance, 60, 0., 300.);
  INT_HISTO(Q2, 50, 0., 10.);
  INT_HISTO(true_contained_length, 101, -10., 1000.);
  INT_HISTO(true_track_multiplicity, 10, 0., 10);
  INT_HISTO(crosses_tpc, 2, -0.5, 1.5);
  INT_HISTO(dist_to_match, 101, -1., 100.);
  INT_HISTO(primary_track_completion, 100, 0., 1.);
  INT_HISTO(n_reco_vertices, 10, -0.5, 9.5);

#undef INT_HISTO
}


void InteractionHistos::Fill(
  unsigned vertex_index,
  bool is_truth,
  const numu::RecoEvent &event,
  const std::vector<event::Interaction> &core_truth) {

  const numu::RecoInteraction &vertex = (is_truth) ? event.truth[vertex_index] : event.reco[vertex_index];

  const std::vector<numu::RecoInteraction> &truth = event.truth;

  const std::map<size_t, numu::RecoTrack> *vertex_tracks;
  if (is_truth) {
    // find the closest reconstructed vertex to this one
    double dist = numu::dist2Match(vertex, event.reco);
    dist_to_match->Fill(dist);
    primary_track_completion->Fill(numu::trackMatchCompletion(vertex_index, event));
    vertex_tracks = &event.true_tracks;
  }
  // closest reconstructed vertex to this one (already contained in object)
  else {
    dist_to_match->Fill(vertex.match.truth_vertex_distance);
    primary_track_completion->Fill(vertex.primary_track.match.completion);
    vertex_tracks = &event.reco_tracks;
  }

  n_reco_vertices->Fill(event.reco.size());

  double track_length_val = vertex.slice.primary_track_index >= 0 ? vertex_tracks->at(vertex.slice.primary_track_index).length: -1;
  track_length->Fill(track_length_val);
  if (vertex.slice.primary_track_index >= 0 && vertex_tracks->at(vertex.slice.primary_track_index).match.has_match) {
    int mcparticle_id = vertex_tracks->at(vertex.slice.primary_track_index).match.mcparticle_id;
  
    double true_track_momentum = event.true_tracks.at(mcparticle_id).momentum; 
    track_p->Fill(true_track_momentum);
    true_deposited_energy->Fill(event.true_tracks.at(mcparticle_id).deposited_energy);

    int crosses_tpc_val = event.true_tracks.at(mcparticle_id).crosses_tpc;
    crosses_tpc->Fill(crosses_tpc_val);

    double length = event.true_tracks.at(mcparticle_id).length;
    true_contained_length->Fill(length);
  }
   
  if (vertex.match.event_track_id > 0 && vertex.primary_track.match.is_primary) {
    int event_id = vertex.match.event_track_id;
    int mctruth_id = vertex.match.mctruth_track_id;

    true_track_multiplicity->Fill(truth[event_id].multiplicity);
    
    if (mctruth_id >= 0) {
      nuE->Fill(core_truth[mctruth_id].neutrino.energy);
      Q2->Fill(core_truth[mctruth_id].neutrino.Q2);
      // get the distance from the beam center
      /*
      float beam_center_distance = sqrt( (core_truth[mctruth_id].neutrino.position.X() - _config.beamCenterX) * 
        (core_truth[mctruth_id].neutrino.position.X() - _config.beamCenterX) +
        (core_truth[mctruth_id].neutrino.position.Y() - _config.beamCenterY) *
        (core_truth[mctruth_id].neutrino.position.Y() - _config.beamCenterY));

      beam_center_distance->Fill(beam_center_distance);
      */
    }
  }
}

  } // namespace SBNOsc
} // namespace ana

