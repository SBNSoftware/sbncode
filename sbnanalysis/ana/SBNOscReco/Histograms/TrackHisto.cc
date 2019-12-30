#include "TrackHisto.h"

#include "TH1D.h"
#include "TH2D.h"

namespace ana {
  namespace SBNOsc {

void TrackHistos::Initialize(const std::string &postfix, const geo::BoxBoundedGeo &detector_volume, double max_length) {
#define TRACK_HISTO(name, n_bins, lo, hi)    name = TH1Shared(new TH1D((#name"_" + postfix).c_str(), #name, n_bins, lo, hi)); fAllHistos.push_back(name.Get())
#define TRACK_2DHISTO(name, binx, lo_x, hi_x, biny, lo_y, hi_y)  name = TH2Shared(new TH2D((#name"_" + postfix).c_str(), #name, binx, lo_x, hi_x, biny, lo_y, hi_y)); fAllHistos.push_back(name.Get())

  TRACK_HISTO(chi2_muon_diff, 100, 0., 1000.);
  TRACK_HISTO(chi2_proton_diff, 100, 0, 1000);
  TRACK_HISTO(chi2_kaon_diff, 100, 0, 1000);
  TRACK_HISTO(chi2_pion_diff, 100, 0, 1000);

  TRACK_HISTO(chi2_muon, 100, 0, 100);
  TRACK_HISTO(chi2_pion, 100, 0, 100);
  TRACK_HISTO(chi2_kaon, 100, 0, 100);
  TRACK_HISTO(chi2_proton, 100, 0, 100);

  TRACK_HISTO(chi2_proton_m_muon, 200, -1000, 1000);

  TRACK_HISTO(range_p, 100, 0., 2.);
  TRACK_HISTO(mcs_p, 100, 0., 2.);
  TRACK_HISTO(deposited_e_max, 100, 0., 2.);
  TRACK_HISTO(deposited_e_avg, 100, 0., 2.);
  
  TRACK_HISTO(range_p_minus_truth, 100, -2., 2);
  TRACK_HISTO(mcs_p_minus_truth, 100, -2., 2.);

  TRACK_2DHISTO(range_p_minus_truth_length, 60, 0., 600., 50, -1., 1);
  TRACK_2DHISTO(mcs_p_minus_truth_length, 60, 0., 600., 50, -1., 1.);

  TRACK_HISTO(deposited_e_max_minus_truth, 100, -2., 2.);
  TRACK_HISTO(deposited_e_avg_minus_truth, 100, -2., 2.);
  TRACK_HISTO(deposited_e_med_minus_truth, 100, -2., 2.); 

  TRACK_HISTO(length, 100, 0., 600.);

  TRACK_HISTO(reco_momentum, 100, 0., 5.);
  TRACK_HISTO(is_contained, 2, -0.5, 1.5);
  
  //TRACK_2DHISTO(range_p_diff, 25, 0, 2.5, 40, -2., 2.); 
  //TRACK_2DHISTO(mcs_p_diff, 25, 0., 2.5, 40, -2., 2.);
  // TRACK_2DHISTO(deposited_e_max_diff, 25, 0., 2.5, 40, -2., 2.);

  TRACK_2DHISTO(range_p_comp, 25, 0, 2.5, 25, 0., 2.5);
  TRACK_2DHISTO(mcs_p_comp, 25, 0., 2.5, 25, 0., 2.5);
  // TRACK_2DHISTO(deposited_e_max_comp,  25, 0., 2.5, 25, 0., 2.5);

  //TRACK_2DHISTO(dQdx_length, 100, 0., 1000., 100, 0., max_length);

  TRACK_HISTO(border_y, 400, detector_volume.MinY(), detector_volume.MaxY());
  TRACK_HISTO(border_x, 400, detector_volume.MinX(), detector_volume.MaxX());
  TRACK_HISTO(border_z, 500, detector_volume.MinZ(), detector_volume.MaxZ()); 
  TRACK_HISTO(true_start_time, 1400, -4000., 3000.);
  TRACK_HISTO(true_start_time_zoom, 3000, -1., 2.);

  TRACK_HISTO(wall_enter, 7, -0.5, 6.5);
  TRACK_HISTO(wall_exit, 7, -0.5, 6.5);
  
  // timing histos
  TRACK_HISTO(has_crt_track_match, 3, -0.5, 1.5);
  TRACK_HISTO(has_crt_hit_match, 3, -0.5, 1.5); 
  TRACK_HISTO(crt_hit_distance, 80, 0., 400.);
  TRACK_HISTO(crt_track_angle, 150, 0., 3.);

  // TRACK_HISTO(flash_match_time, 2000, -0.2, 1.8);
  // TRACK_HISTO(crt_v_flash_match_time, 2000, -4., 4.);
  
  double min_matchtime_t = -1640;
  double max_matchtime_t =  3280;
  int n_matchtime_bins = 1000;
  
  double min_comptime = -0.5;
  double max_comptime = 0.5;
  int n_comptime_bins = 1000;
  
  TRACK_HISTO(crt_match_time, n_matchtime_bins, min_matchtime_t, max_matchtime_t);

  TRACK_HISTO(completion, 200, -1, 1);

  TRACK_HISTO(stopping_chisq_start, 100, 0., 10.);
  TRACK_HISTO(stopping_chisq_finish, 100, 0., 10.);
  TRACK_HISTO(stopping_chisq, 100., 0., 10.);

  TRACK_2DHISTO(pid_confusion_tr, 2, -0.5, 1.5, 2, -0.5, 1.5);

#undef TRACK_HISTO
#undef TRACK_2DHISTO
}

void TrackHistos::Fill(
    const numu::RecoTrack &track, 
    const std::map<size_t, numu::RecoTrack> &true_tracks) {
#define FILL(hist, x) hist.Fill(x)
#define FILL2D(hist, x, y) hist.Fill(x, y)

  // Primary track histos
  if (track.min_chi2 > 0) {
    FILL(chi2_proton_diff, track.chi2_proton - track.min_chi2);
    FILL(chi2_muon_diff, track.chi2_muon - track.min_chi2);
    FILL(chi2_pion_diff, track.chi2_pion - track.min_chi2);
    FILL(chi2_kaon_diff, track.chi2_kaon - track.min_chi2);

    FILL(chi2_muon, track.chi2_muon/track.pid_n_dof);
    FILL(chi2_kaon, track.chi2_kaon/track.pid_n_dof);
    FILL(chi2_pion, track.chi2_pion/track.pid_n_dof);
    FILL(chi2_proton, track.chi2_proton/track.pid_n_dof);
    FILL(chi2_proton_m_muon, track.chi2_proton - track.chi2_muon);

    bool is_proton_reco = track.chi2_proton < track.chi2_muon;
    if (track.match.has_match) {
      bool is_proton_true = abs(track.match.match_pdg) == 2212;
      bool is_muon_true = abs(track.match.match_pdg) == 13;
      if (is_proton_true || is_muon_true) {
        FILL2D(pid_confusion_tr, is_proton_true, is_proton_reco);
      }
    }
  }
	  
  FILL(range_p, track.range_momentum);
  FILL(mcs_p, track.mcs_momentum);
  FILL(deposited_e_max, track.deposited_energy_max);
  
  FILL(length, track.length);

  FILL(reco_momentum, track.momentum);
  FILL(is_contained, track.is_contained);

  //dQdx_length->Fill(track.mean_trucated_dQdx, track.length);

  if (std::min(abs(track.start.Y() - border_y.Get()->GetBinLowEdge(1)), 
               abs(track.start.Y() - (border_y.Get()->GetBinLowEdge(border_y.Get()->GetNbinsX()) + border_y.Get()->GetBinWidth(border_y.Get()->GetNbinsX()))))
   <  std::min(abs(track.end.Y() - border_y.Get()->GetBinLowEdge(1)), 
               abs(track.end.Y() - (border_y.Get()->GetBinLowEdge(border_y.Get()->GetNbinsX()) + border_y.Get()->GetBinWidth(border_y.Get()->GetNbinsX()))))) {
    FILL(border_y, track.start.Y());
  }
  else {
    FILL(border_y, track.end.Y());
  }

  if (std::min(abs(track.start.X() - border_x.Get()->GetBinLowEdge(1)), 
               abs(track.start.X() - (border_x.Get()->GetBinLowEdge(border_x.Get()->GetNbinsX()) + border_x.Get()->GetBinWidth(border_x.Get()->GetNbinsX()))))
   <  std::min(abs(track.end.X() - border_x.Get()->GetBinLowEdge(1)), 
               abs(track.end.X() - (border_x.Get()->GetBinLowEdge(border_x.Get()->GetNbinsX()) + border_x.Get()->GetBinWidth(border_x.Get()->GetNbinsX()))))) {
    FILL(border_x, track.start.X());
  }
  else {
    FILL(border_x, track.end.X());
  }

  if (std::min(abs(track.start.Z() - border_z.Get()->GetBinLowEdge(1)), 
               abs(track.start.Z() - (border_z.Get()->GetBinLowEdge(border_z.Get()->GetNbinsX()) + border_z.Get()->GetBinWidth(border_z.Get()->GetNbinsX()))))
   <  std::min(abs(track.end.Z() - border_z.Get()->GetBinLowEdge(1)), 
               abs(track.end.Z() - (border_z.Get()->GetBinLowEdge(border_z.Get()->GetNbinsX()) + border_z.Get()->GetBinWidth(border_z.Get()->GetNbinsX()))))) {
    FILL(border_z, track.start.Z());
  }
  else {
    FILL(border_z, track.end.Z());
  }

  FILL(has_crt_hit_match, track.crt_match.hit_match.present);
  FILL(has_crt_track_match, track.crt_match.track.present);
  if (track.crt_match.track.present) {
    FILL(crt_match_time, track.crt_match.track.time);
    FILL(crt_track_angle, track.crt_match.track.angle);
  }
  else if (track.crt_match.hit_match.present) {
    FILL(crt_match_time, track.crt_match.hit_match.time);
    FILL(crt_hit_distance, track.crt_match.hit_match.distance);
  }

  /*
  has_flash_match->Fill(track.flash_match.present);
  if (track.flash_match.present) {
    const numu::FlashMatch &flash_match = track.flash_match;
    double flash_time = flash_match.match_time_first;
    flash_match_time->Fill(flash_time);
    if (track.crt_match.track.present) {
      crt_v_flash_match_time->Fill(track.crt_match.track.time - flash_time);
    }
    else if (track.crt_match.hit_match.present) {
      crt_v_flash_match_time->Fill(track.crt_match.hit_match.time - flash_time);
    }  
  }*/
  
  // check if truth match
  if (track.match.has_match && track.match.mcparticle_id >= 0) {
    const numu::RecoTrack &true_track = true_tracks.at(track.match.mcparticle_id);
    FILL(range_p_minus_truth, (track.range_momentum - true_track.momentum) / true_track.momentum);
    FILL(mcs_p_minus_truth, (track.mcs_momentum - true_track.momentum) / true_track.momentum);

    FILL2D(range_p_minus_truth_length, track.length, (track.range_momentum - true_track.momentum) / true_track.momentum);
    FILL2D(mcs_p_minus_truth_length, track.length, (track.mcs_momentum - true_track.momentum) / true_track.momentum);

    FILL(deposited_e_max_minus_truth, track.deposited_energy_max - true_track.energy);
    FILL(deposited_e_avg_minus_truth, track.deposited_energy_avg - true_track.energy);
    FILL(deposited_e_med_minus_truth, track.deposited_energy_med - true_track.energy);
    
    //range_p_diff->Fill(true_track.momentum, track.range_momentum - true_track.momentum);
    //mcs_p_diff->Fill(true_track.momentum, track.mcs_momentum - true_track.momentum);
    // deposited_e_max_diff->Fill(true_track.energy, track.deposited_energy_max - true_track.energy);
    
    FILL2D(range_p_comp, true_track.momentum, track.range_momentum);
    FILL2D(mcs_p_comp, true_track.momentum, track.mcs_momentum);
    // deposited_e_max_comp->Fill(true_track.energy, track.deposited_energy_max);

    FILL(completion, track.match.completion);

    FILL(wall_enter, true_track.wall_enter);
    FILL(wall_exit, true_track.wall_exit);
    FILL(true_start_time, true_track.start_time);
    FILL(true_start_time_zoom, true_track.start_time);
  }
  else {
    FILL(completion, -0.5);
  }

  FILL(stopping_chisq_start, track.stopping_chisq_start);
  FILL(stopping_chisq_finish, track.stopping_chisq_finish);

  FILL(stopping_chisq, track.stopping_chisq_start);
  FILL(stopping_chisq, track.stopping_chisq_finish);

#undef FILL
#undef FILL2D
}
  } // namespace SBNOsc
} // namespace ana
