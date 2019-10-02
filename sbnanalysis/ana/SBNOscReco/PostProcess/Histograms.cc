#include <TH1D.h>
#include <TH2D.h>

#include "nusimdata/SimulationBase/MCTruth.h"

#include "Histograms.h"

namespace ana {
  namespace SBNOsc {

void FillTrack(
  const numu::RecoTrack &track, 
  const std::map<size_t, numu::RecoTrack> &true_tracks, 
  TrackHistos* track_histos) {

  std::vector<bool> do_fill {true, false, false, false};
  if (track.match.has_match) {
    do_fill[1] = track.match.mctruth_origin == simb::Origin_t::kCosmicRay;
    do_fill[2] = track.match.mctruth_origin == simb::Origin_t::kBeamNeutrino;
  }
  if (!do_fill[1] && !do_fill[2]) do_fill[3] = true;
  for (int i = 0; i < do_fill.size(); i++) {
    if (do_fill[i]) {
      track_histos[i].Fill(track, true_tracks);
    }
  }
}

void Histograms::Fill(const numu::RecoEvent &event, const event::Event &core, const std::array<bool, Cuts::nCuts> &cuts) {
  for (const numu::RecoTrack &track: event.tracks) {
    FillTrack(track, event.true_tracks, fAllTracks);
  }

  for (const numu::RecoInteraction &interaction: event.reco) {
    if (event.tracks.size() > (unsigned)interaction.slice.primary_track_index) {
      FillTrack(event.tracks.at(interaction.slice.primary_track_index), event.true_tracks, fPrimaryTracks);
    }
    // fill histos
    for (size_t cut_i=0; cut_i < Cuts::nCuts; cut_i++) {
      int mode = interaction.match.mode; 
      if (cuts[cut_i]) {
        fInteraction[cut_i+InteractionHistos::recoCutOffset][mode].Fill(interaction, event.truth, core.truth);
        fInteraction[cut_i+InteractionHistos::recoCutOffset][numu::mAll].Fill(interaction, event.truth, core.truth);
      }
    }
  }

  for (const numu::RecoInteraction &truth: event.truth) {
    // TODO: fix
    std::array<bool, 2> fill_truth = {true, event.reco.size() >= 1};
    for (int cut_i = 0; cut_i < 2; cut_i++) {
      int mode = truth.match.mode;
      fInteraction[cut_i][mode].Fill(truth, event.truth, core.truth);
      fInteraction[cut_i][numu::mAll].Fill(truth, event.truth, core.truth);
    }
  }
}

void InteractionHistos::Write() {
  track_length->Write();
  track_p->Write();
  nuE->Write();
  beam_center_distance->Write();
  Q2->Write();
  true_contained_length->Write();
  true_track_multiplicity->Write();
  crosses_tpc->Write();
}

void TrackHistos::Write() {
  chi2_proton_diff->Write();
  chi2_muon_diff->Write();
  chi2_pion_diff->Write();
  chi2_kaon_diff->Write();
  
  range_p->Write();
  mcs_p->Write();
  deposited_e_max->Write();
  deposited_e_avg->Write();
  
  range_p_minus_truth->Write();
  mcs_p_minus_truth->Write();
  deposited_e_max_minus_truth->Write();
  deposited_e_avg_minus_truth->Write();
  deposited_e_med_minus_truth->Write();
  
  length->Write();
  is_contained->Write();
  
  range_p_diff->Write();
  mcs_p_diff->Write();
  deposited_e_max_diff->Write();
  
  range_p_comp->Write();
  mcs_p_comp->Write();
  deposited_e_max_comp->Write();
  
  has_crt_track_match->Write();
  has_crt_hit_match->Write();
  has_flash_match->Write();
  crt_hit_match_time->Write();
  flash_match_time->Write();
  crt_v_flash_match_time->Write();
}

Histograms::Histograms() {
  for (unsigned i = 0; i < InteractionHistos::nHistos; i++) {
    for (const auto mode: InteractionHistos::allModes) {
      fInteraction[i][mode].Initialize("", mode, i); 
    }
  }

  for (unsigned i = 0; i < TrackHistos::nTrackHistos; i++) {
    fAllTracks[i].Initialize("All", i);
    fPrimaryTracks[i].Initialize("Primary", i);
  } 
}

void Histograms::Write() {
  for (unsigned i = 0; i < InteractionHistos::nHistos; i++) {
    for (const auto mode: InteractionHistos::allModes) {
      fInteraction[i][mode].Write();
    }
  }

  for (unsigned i = 0; i < TrackHistos::nTrackHistos; i++) {
    fAllTracks[i].Write();
    fPrimaryTracks[i].Write();
  } 
}

void InteractionHistos::Initialize(const std::string &prefix, numu::InteractionMode mode, unsigned i) {
  track_length = new TH1D(("track_length_" + mode2Str(mode) + "_" + prefix + histoNames[i]).c_str(), "track_length", 101, -10, 1000);
  track_p = new TH1D(("track_p_" + mode2Str(mode) + "_" + prefix + histoNames[i]).c_str(), "track_p", 50, 0., 5.);
  nuE = new TH1D(("nuE_" + mode2Str(mode) +"_" + prefix + histoNames[i]).c_str(), "nuE", 50, 0., 5.);
  beam_center_distance = new TH1D(("beam_dist_" + mode2Str(mode) + "_" + prefix + histoNames[i]).c_str(), "beam_dist", 60, 0., 300.);
  Q2 = new TH1D(("Q2_" + mode2Str(mode) + "_" + prefix + histoNames[i]).c_str(), "Q2", 50, 0., 10.);
  true_contained_length = new TH1D(("tt_contained_length_" + mode2Str(mode) + "_" + prefix + histoNames[i]).c_str(), "tt_contained_length", 101, -10., 1000.);
  true_track_multiplicity = new TH1D(("true_track_multiplicity_" + mode2Str(mode) + "_" + prefix + histoNames[i]).c_str(), "true_track_multiplicity", 10, 0., 10.);
  crosses_tpc = new TH1D(("crosses_tpc_" + mode2Str(mode) + "_" + prefix + histoNames[i]).c_str(), "crosses_tpc", 2, -0.5, 1.5);
}

void TrackHistos::Initialize(const std::string &prefix, unsigned i) {
  chi2_muon_diff = new TH1D((std::string("chi2_muon_diff_") + prefix + trackHistoNames[i]).c_str(), "chi2_muon_diff", 100, 0., 100.);
  
  chi2_proton_diff = new TH1D((std::string("chi2_proton_diff_") + prefix + trackHistoNames[i]).c_str(), "chi2_proton_diff", 101, -0.1, 10);
  chi2_kaon_diff = new TH1D((std::string("chi2_kaon_diff_") + prefix + trackHistoNames[i]).c_str(), "chi2_kaon_diff", 101, -0.1, 10);
  chi2_pion_diff = new TH1D((std::string("chi2_pion_diff_") + prefix + trackHistoNames[i]).c_str(), "chi2_pion_diff", 101, -0.1, 10);
  
  range_p = new TH1D((std::string("range_p_") + prefix + trackHistoNames[i]).c_str(), "range_p", 100, 0., 2.);
  mcs_p = new TH1D((std::string("mcs_p_") + prefix + trackHistoNames[i]).c_str(), "mcs_p", 100, 0., 2.);
  deposited_e_max = new TH1D((std::string("deposited_e_max_") + prefix + trackHistoNames[i]).c_str(), "deposited_e_max", 100, 0., 2.);
  deposited_e_avg = new TH1D((std::string("deposited_e_avg_") + prefix + trackHistoNames[i]).c_str(), "deposited_e_avg", 100, 0., 2.);
  
  range_p_minus_truth = new TH1D((std::string("range_p_minus_truth_") + prefix + trackHistoNames[i]).c_str(), "range_p_minus_truth", 100, -2., 2.);
  mcs_p_minus_truth = new TH1D((std::string("mcs_p_minus_truth_") + prefix + trackHistoNames[i]).c_str(), "mcs_p_minus_truth", 100, -2., 2.);
  deposited_e_max_minus_truth = new TH1D((std::string("deposited_e_max_minus_truth_") + prefix + trackHistoNames[i]).c_str(), "deposited_e_max_minus_truth", 100, -2., 2.);
  deposited_e_avg_minus_truth = new TH1D((std::string("deposited_e_avg_minus_truth_") + prefix + trackHistoNames[i]).c_str(), "deposited_e_avg_minus_truth", 100, -2., 2.);
  deposited_e_med_minus_truth = new TH1D((std::string("deposited_e_med_minus_truth_") + prefix + trackHistoNames[i]).c_str(), "deposited_e_med_minus_truth", 100, -2., 2.);
  
  length = new TH1D((std::string("length_") + prefix + trackHistoNames[i]).c_str(), "length", 100, 0., 500.);
  is_contained = new TH1D((std::string("is_contained_") + prefix + trackHistoNames[i]).c_str(), "is_contained", 2, -0.5, 1.5);
  
  
  range_p_diff = new TH2D((std::string("range_p_diff_") + prefix + trackHistoNames[i]).c_str(), "range_p_diff", 25, 0, 2.5, 40, -2., 2.);
  mcs_p_diff = new TH2D((std::string("mcs_p_diff_") + prefix + trackHistoNames[i]).c_str(), "mcs_p_diff", 25, 0., 2.5, 40, -2., 2.);
  deposited_e_max_diff = new TH2D((std::string("deposited_e_max_diff_") + prefix + trackHistoNames[i]).c_str(), "deposited_e_max_diff", 25, 0., 2.5, 40, -2., 2.);
  
  
  range_p_comp = new TH2D((std::string("range_p_comp_") + prefix + trackHistoNames[i]).c_str(), "range_p_comp", 25, 0, 2.5, 25, 0., 2.5);
  mcs_p_comp = new TH2D((std::string("mcs_p_comp_") + prefix + trackHistoNames[i]).c_str(), "mcs_p_comp", 25, 0., 2.5, 25, 0., 2.5);
  deposited_e_max_comp = new TH2D((std::string("deposited_e_max_comp_") + prefix + trackHistoNames[i]).c_str(), "deposited_e_max_comp", 25, 0., 2.5, 25, 0., 2.5);
  
  // timing histos
  has_crt_track_match = new TH1D((std::string("has_crt_track_match_") + prefix + trackHistoNames[i]).c_str(), "has_crt_track_match", 3, -0.5, 1.5);
  has_crt_hit_match = new TH1D((std::string("has_crt_hit_match_") + prefix + trackHistoNames[i]).c_str(), "has_crt_hit_match", 3, -0.5, 1.5);
  has_flash_match = new TH1D((std::string("has_flash_match_") + prefix + trackHistoNames[i]).c_str(), "has_flash_match", 3, -0.5, 1.5);
  
  double min_matchtime_t = -1640;
  double max_matchtime_t =  3280;
  int n_matchtime_bins = 1000;
  
  double min_comptime = -0.5;
  double max_comptime = 0.5;
  int n_comptime_bins = 1000;
  
  crt_hit_match_time = new TH1D((std::string("crt_hit_match_time_") + prefix + trackHistoNames[i]).c_str(), "crt_hit_match_time", n_matchtime_bins, min_matchtime_t, max_matchtime_t);
  flash_match_time = new TH1D((std::string("flash_match_time_") + prefix + trackHistoNames[i]).c_str(), "flash_match_time", n_matchtime_bins, min_matchtime_t, max_matchtime_t);
  crt_v_flash_match_time = new TH1D((std::string("crt_v_flash_match_time_") + prefix + trackHistoNames[i]).c_str(), "crt_v_flash_match_time", n_comptime_bins, min_comptime, max_comptime);
}

void TrackHistos::Fill(
    const numu::RecoTrack &track, 
    const std::map<size_t, numu::RecoTrack> &true_tracks) {

  // Primary track histos
  if (track.min_chi2 > 0) {
    chi2_proton_diff->Fill(track.chi2_proton - track.min_chi2);
    chi2_muon_diff->Fill(track.chi2_muon - track.min_chi2);
    chi2_pion_diff->Fill(track.chi2_pion - track.min_chi2);
    chi2_kaon_diff->Fill(track.chi2_kaon - track.min_chi2);
  }
	  
  range_p->Fill(track.range_momentum_muon); 
  double mcs_p_val = track.mcs_is_backward ? track.bwd_mcs_momentum_muon : track.fwd_mcs_momentum_muon;
  mcs_p->Fill(mcs_p_val);
  deposited_e_max->Fill(track.deposited_energy_max);
  
  length->Fill(track.length);
  is_contained->Fill(track.is_contained);
  
  bool has_crt_match = track.crt_match.size() > 0;
  bool has_flash_match_val = track.flash_match.size() > 0;
  has_flash_match->Fill(has_flash_match_val);
  double crt_match_time;
  bool has_crt_match_time = false;
  if (has_crt_match) {
    const numu::CRTMatch &crt_match = track.crt_match.at(0);
    if (crt_match.has_track_match) {
      has_crt_track_match->Fill(1.);
      has_crt_hit_match->Fill(1.);
      crt_match_time = crt_match.match_time;
      has_crt_match_time = true;
    }
    else if (crt_match.has_hit_match) {
      has_crt_track_match->Fill(0.);
      has_crt_hit_match->Fill(1.);
      crt_match_time = crt_match.match_time;
      has_crt_match_time = true;
      crt_hit_match_time->Fill(crt_match_time);
    }
    else {
      has_crt_track_match->Fill(0.);
      has_crt_hit_match->Fill(0.);
    }
  }
  else {
    has_crt_track_match->Fill(0.);
    has_crt_hit_match->Fill(0.);
  }
  if (has_flash_match_val) {
    const numu::FlashMatch &flash_match = track.flash_match.at(0);
    double flash_time = flash_match.match_time_first;
    std::cout << "FLASH match time: " << flash_time << std::endl;
    flash_match_time->Fill(flash_time);
    if (has_crt_match_time) {
      crt_v_flash_match_time->Fill(crt_match_time - flash_time);  
      std::cout << "CRT match time: " << crt_match_time << std::endl;
    }
  }
  
  // check if truth match
  if (track.match.has_match && track.match.mcparticle_id >= 0) {
    const numu::RecoTrack &true_track = true_tracks.at(track.match.mcparticle_id);
    range_p_minus_truth->Fill(track.range_momentum_muon - true_track.momentum);
    mcs_p_minus_truth->Fill(mcs_p_val - true_track.momentum); 
    deposited_e_max_minus_truth->Fill(track.deposited_energy_max - true_track.energy);
    deposited_e_avg_minus_truth->Fill(track.deposited_energy_avg - true_track.energy);
    deposited_e_med_minus_truth->Fill(track.deposited_energy_med - true_track.energy);
    
    //std::cout << "AT FILL -- is contained: " << true_track.is_contained << std::endl;      
    
    range_p_diff->Fill(true_track.momentum, track.range_momentum_muon - true_track.momentum);
    mcs_p_diff->Fill(true_track.momentum, mcs_p_val - true_track.momentum);
    deposited_e_max_diff->Fill(true_track.energy, track.deposited_energy_max - true_track.energy);
    
    range_p_comp->Fill(true_track.momentum, track.range_momentum_muon);
    mcs_p_comp->Fill(true_track.momentum, mcs_p_val);
    deposited_e_max_comp->Fill(true_track.energy, track.deposited_energy_max);
  }
}

void InteractionHistos::Fill(
  const numu::RecoInteraction &vertex, 
  const std::vector<numu::RecoInteraction> &truth, 
  const std::vector<event::Interaction> &core_truth) {

  double track_length_val = vertex.slice.primary_track_index >= 0 ? vertex.slice.tracks.at(vertex.slice.primary_track_index).length: -1;
  track_length->Fill(track_length_val);
  if (vertex.match.event_track_id >= 0) {
    int event_id = vertex.match.event_track_id;
    int mctruth_id = vertex.match.mctruth_track_id;
    
    double true_track_momentum = truth[event_id].slice.primary_track_index >= 0 ? truth[event_id].slice.tracks.at(truth[event_id].slice.primary_track_index).momentum : -1;
    track_p->Fill(true_track_momentum);

    int crosses_tpc_val = truth[event_id].slice.primary_track_index >= 0 ? truth[event_id].slice.tracks.at(truth[event_id].slice.primary_track_index).crosses_tpc: -1;
    crosses_tpc->Fill(crosses_tpc_val);

    double length = truth[event_id].slice.primary_track_index >= 0 ? truth[event_id].slice.tracks.at(truth[event_id].slice.primary_track_index).length: -1;
    true_contained_length->Fill(length);
    
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
