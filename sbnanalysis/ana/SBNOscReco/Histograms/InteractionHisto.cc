#include "InteractionHisto.h"
#include "Derived.h"
#include "../TriggerEmulator/PMTTrigger.h"

#include "TH1D.h"
#include "TH2D.h"

namespace ana {
  namespace SBNOsc {

void InteractionHistos::Initialize(const std::string &postfix, const geo::BoxBoundedGeo &detector_volume, const std::vector<double> &tagger_volume) {
#define INT_HISTO(name, n_bins, lo, hi)    name = TH1Shared(new TH1D((#name"_" + postfix).c_str(), #name, n_bins, lo, hi)); fAllHistos.push_back(name.Get())
#define INT_HISTO2D(name, n_binsx, xlo, xhi, n_binsy, ylo, yhi) name = TH2Shared(new TH2D((#name"_" + postfix).c_str(), #name, n_binsx, xlo, xhi, n_binsy, ylo, yhi)); fAllHistos.push_back(name.Get())
#define INT_HISTO2D_BINSY(name, n_binsx, xlo, xhi, n_binsy, binsy) name = TH2Shared(new TH2D((#name"_" + postfix).c_str(), #name, n_binsx, xlo, xhi, n_binsy, binsy)); fAllHistos.push_back(name.Get())

  INT_HISTO(track_length, 100, 0., 600.);
  INT_HISTO(track_p, 50, 0., 5.);
  INT_HISTO(true_deposited_energy, 50., 0., 5.);
  INT_HISTO(nuE, 50, 0., 5.);
  INT_HISTO(beam_center_distance, 60, 0., 300.);
  INT_HISTO(Q2, 50, 0., 10.);
  INT_HISTO(true_contained_length, 100, 0., 600.);
  INT_HISTO(true_track_multiplicity, 10, 0., 10);
  INT_HISTO(crosses_tpc, 2, -0.5, 1.5);
  INT_HISTO(dist_to_match, 101, -1., 100.);
  INT_HISTO(primary_track_completion, 100, 0., 1.);
  INT_HISTO(n_reco_vertices, 10, -0.5, 9.5);
  INT_HISTO(maxpe_crt_intime_hit, 1000, 0., 1000.);
  INT_HISTO(crt_pes, 1000, 0., 1000.);

  INT_HISTO(crt_hit_times, 300, -20., 10.);
  INT_HISTO(closest_crt_hit_time, 300, -20., 10.);

  INT_HISTO(fmatch_score, 200, 0., 200.);
  INT_HISTO2D(fmatch_score_true_time, 100, 0., 100., 1400, -4000., 3000.);
  INT_HISTO2D(fmatch_score_true_time_zoom, 100, 0., 100., 100, -5., 5.);
  double intout_times[4] = {-3000., 0., 1.6, 3000.};
  INT_HISTO(fmatch_score_intime, 100, 0., 100.);
  INT_HISTO(fmatch_score_outtime, 100, 0., 100.);
  INT_HISTO2D(fmatch_time_true_time_zoom, 100, -5., 5., 100, -5., 5.);
  INT_HISTO(fmatch_time_real_time, 4000, -2, 2);
  INT_HISTO(fmatch_time, 700, -5., 2.);

  // INT_HISTO2D(light_trigger, 20, 0.5, 20.5, 200, 6000, 8000);

  INT_HISTO2D(intime_crt_hits_xy, 100, tagger_volume[0], tagger_volume[3], 100, tagger_volume[1], tagger_volume[4]);
  INT_HISTO2D(intime_crt_hits_xz, 100, tagger_volume[0], tagger_volume[3], 100, tagger_volume[2], tagger_volume[5]);
  INT_HISTO2D(intime_crt_hits_yz, 100, tagger_volume[1], tagger_volume[4], 100, tagger_volume[2], tagger_volume[5]);

  INT_HISTO2D(vertex_xy, 100, detector_volume.MinX(), detector_volume.MaxX(), 100, detector_volume.MinY(), detector_volume.MaxY());
  INT_HISTO2D(vertex_yz, 100, detector_volume.MinY(), detector_volume.MaxY(), 100, detector_volume.MinZ(), detector_volume.MaxZ());
  INT_HISTO2D(vertex_xz, 100, detector_volume.MinX(), detector_volume.MaxX(), 100, detector_volume.MinZ(), detector_volume.MaxZ());

#undef INT_HISTO
}


void InteractionHistos::Fill(
  unsigned vertex_index,
  bool is_truth,
  const numu::RecoEvent &event,
  const std::vector<event::Interaction> &core_truth) {
#define FILL(hist, val) hist.Fill(val);
#define FILL2D(hist, x, y) hist.Fill(x, y);

  const numu::RecoInteraction &vertex = (is_truth) ? event.truth[vertex_index] : event.reco[vertex_index];

  const std::vector<numu::RecoInteraction> &truth = event.truth;

  const std::map<size_t, numu::RecoTrack> *vertex_tracks;
  if (is_truth) {
    // find the closest reconstructed vertex to this one
    double dist = numu::dist2Match(vertex, event.reco);
    FILL(dist_to_match, dist);
    FILL(primary_track_completion, numu::trackMatchCompletion(vertex_index, event));
    vertex_tracks = &event.true_tracks;
  }
  // closest reconstructed vertex to this one (already contained in object)
  else {
    FILL(dist_to_match, vertex.match.truth_vertex_distance);
    FILL(primary_track_completion, vertex.primary_track.match.completion);
    vertex_tracks = &event.reco_tracks;
  }

  //std::vector<int> thresholds = numu::TriggerThresholds(event.flash_trigger_primitives, 20);
  //for (int i = 0; i < thresholds.size(); i++) {
  //  light_trigger->Fill(i, thresholds[i]);
  //}

  FILL2D(vertex_xy, vertex.position.X(), vertex.position.Y());
  FILL2D(vertex_yz, vertex.position.Y(), vertex.position.Z());
  FILL2D(vertex_xz, vertex.position.X(), vertex.position.Z());

  FILL(n_reco_vertices, event.reco.size());

  double maxpe = 0.;
  double closest_time_dist = -1;
  double closest_time = 0.;
  for (const numu::CRTHit &hit: event.in_time_crt_hits) {
    FILL2D(intime_crt_hits_xy, hit.location.X(), hit.location.Y());
    FILL2D(intime_crt_hits_xz, hit.location.X(), hit.location.Z());
    FILL2D(intime_crt_hits_yz, hit.location.Y(), hit.location.Z());

    FILL(crt_pes, hit.pes);
    if (hit.pes > maxpe) maxpe = hit.pes;


    //if (hit.pes < 100.) continue;

    FILL(crt_hit_times, hit.time);

    if (closest_time_dist < 0. || closest_time_dist > 1e-3) {
       double this_time_dist = -1;
       if (hit.time > -0.2 && hit.time < 1.8) this_time_dist = 0.;
       else if (hit.time < 0.) this_time_dist = -hit.time -0.2;
       if (this_time_dist >= 0. && (this_time_dist < closest_time_dist || closest_time_dist < 0.)) {
         closest_time = hit.time;
         closest_time_dist = this_time_dist;
       }
    }
  }
  if (closest_time_dist >= 0.) {
    FILL(closest_crt_hit_time, closest_time);
  }
  else {
    FILL(closest_crt_hit_time, 30.);
  }

  std::cout << "Max PE: " << maxpe << std::endl; 
  if (event.in_time_crt_hits.size() == 0) {
    FILL(maxpe_crt_intime_hit, -1);
  }
  else {
    FILL(maxpe_crt_intime_hit, maxpe);
  }

  if (vertex.slice.flash_match.present) {
    FILL(fmatch_score, vertex.slice.flash_match.score);
    FILL(fmatch_time, vertex.slice.flash_match.time);
    if (vertex.slice.primary_track_index >= 0 && vertex_tracks->at(vertex.slice.primary_track_index).match.has_match) {
      int mcparticle_id = vertex_tracks->at(vertex.slice.primary_track_index).match.mcparticle_id;
      double true_time = event.true_tracks.at(mcparticle_id).start_time;
      FILL2D(fmatch_score_true_time, vertex.slice.flash_match.score, true_time);
      FILL2D(fmatch_score_true_time_zoom, vertex.slice.flash_match.score, true_time);
      FILL2D(fmatch_time_true_time_zoom, vertex.slice.flash_match.time, true_time);

      if (true_time < 0. || true_time > 1.6) {
        FILL(fmatch_score_outtime, vertex.slice.flash_match.score); 
      }
      else {
        FILL(fmatch_score_intime, vertex.slice.flash_match.score); 
      }

      FILL(fmatch_time_real_time, vertex.slice.flash_match.time - true_time);
    }
  }


  double track_length_val = vertex.slice.primary_track_index >= 0 ? vertex_tracks->at(vertex.slice.primary_track_index).length: -1;
  FILL(track_length, track_length_val);
  if (vertex.slice.primary_track_index >= 0 && vertex_tracks->at(vertex.slice.primary_track_index).match.has_match) {
    int mcparticle_id = vertex_tracks->at(vertex.slice.primary_track_index).match.mcparticle_id;
  
    double true_track_momentum = event.true_tracks.at(mcparticle_id).momentum; 
    FILL(track_p, true_track_momentum);
    FILL(true_deposited_energy, event.true_tracks.at(mcparticle_id).deposited_energy);

    int crosses_tpc_val = event.true_tracks.at(mcparticle_id).crosses_tpc;
    FILL(crosses_tpc, crosses_tpc_val);

    double length = event.true_tracks.at(mcparticle_id).length;
    FILL(true_contained_length, length);
  }
   
  if (vertex.match.event_track_id >= 0) {
    int event_id = vertex.match.event_track_id;
    int mctruth_id = vertex.match.mctruth_track_id;

    FILL(true_track_multiplicity, truth[event_id].multiplicity);
    
    if (mctruth_id >= 0) {
      FILL(nuE, core_truth[mctruth_id].neutrino.energy);
      FILL(Q2, core_truth[mctruth_id].neutrino.Q2);
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
#undef FILL
#undef FILL2D
}

  } // namespace SBNOsc
} // namespace ana

