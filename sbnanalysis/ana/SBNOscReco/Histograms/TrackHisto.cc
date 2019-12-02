#include "TrackHisto.h"

#include "TH1D.h"
#include "TH2D.h"

namespace ana {
  namespace SBNOsc {

void TrackHistos::Initialize(const std::string &postfix, const geo::BoxBoundedGeo &detector_volume, double max_length) {
#define TRACK_HISTO(name, n_bins, lo, hi)    name = new TH1D((#name"_" + postfix).c_str(), #name, n_bins, lo, hi); fAllHistos.push_back(name)
#define TRACK_2DHISTO(name, binx, lo_x, hi_x, biny, lo_y, hi_y)  name = new TH2D((#name"_" + postfix).c_str(), #name, binx, lo_x, hi_x, biny, lo_y, hi_y); fAllHistos.push_back(name)

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
  TRACK_HISTO(deposited_e_max_minus_truth, 100, -2., 2.);
  TRACK_HISTO(deposited_e_avg_minus_truth, 100, -2., 2.);
  TRACK_HISTO(deposited_e_med_minus_truth, 100, -2., 2.); 

  TRACK_HISTO(length, 100, 0., max_length);

  TRACK_HISTO(reco_momentum, 100, 0., 5.);
  TRACK_HISTO(is_contained, 2, -0.5, 1.5);
  
  TRACK_2DHISTO(range_p_diff, 25, 0, 2.5, 40, -2., 2.); 
  TRACK_2DHISTO(mcs_p_diff, 25, 0., 2.5, 40, -2., 2.);
  TRACK_2DHISTO(deposited_e_max_diff, 25, 0., 2.5, 40, -2., 2.);

  TRACK_2DHISTO(range_p_comp, 25, 0, 2.5, 25, 0., 2.5);
  TRACK_2DHISTO(mcs_p_comp, 25, 0., 2.5, 25, 0., 2.5);
  TRACK_2DHISTO(deposited_e_max_comp,  25, 0., 2.5, 25, 0., 2.5);

  TRACK_2DHISTO(dQdx_length, 100, 0., 1000., 100, 0., max_length);

  TRACK_HISTO(border_y, 400, detector_volume.MinY(), detector_volume.MaxY());
  TRACK_HISTO(border_x, 400, detector_volume.MinX(), detector_volume.MaxX());
  TRACK_HISTO(border_z, 500, detector_volume.MinZ(), detector_volume.MaxZ()); 
  TRACK_HISTO(true_start_time, 1400, -4000., 3000.);

  TRACK_HISTO(wall_enter, 7, -0.5, 6.5);
  TRACK_HISTO(wall_exit, 7, -0.5, 6.5);
  
  // timing histos
  TRACK_HISTO(has_crt_track_match, 3, -0.5, 1.5);
  TRACK_HISTO(has_crt_hit_match, 3, -0.5, 1.5); 
  TRACK_HISTO(has_flash_match, 3, -0.5, 1.5);
  TRACK_HISTO(crt_hit_distance, 500, 0., 2000.);
  TRACK_HISTO(crt_track_angle, 150, 0., 3.);
  
  double min_matchtime_t = -1640;
  double max_matchtime_t =  3280;
  int n_matchtime_bins = 1000;
  
  double min_comptime = -0.5;
  double max_comptime = 0.5;
  int n_comptime_bins = 1000;
  
  TRACK_HISTO(crt_match_time, n_matchtime_bins, min_matchtime_t, max_matchtime_t);
  TRACK_HISTO(flash_match_time, n_matchtime_bins, min_matchtime_t, max_matchtime_t);
  TRACK_HISTO(crt_v_flash_match_time, n_comptime_bins, min_comptime, max_comptime);

  TRACK_HISTO(completion, 200, -1, 1);

  TRACK_HISTO(stopping_chisq_start, 100, 0., 10.);
  TRACK_HISTO(stopping_chisq_finish, 100, 0., 10.);
  TRACK_HISTO(stopping_chisq, 100., 0., 10.);

#undef TRACK_HISTO
#undef TRACK_2DHISTO
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

    chi2_muon->Fill(track.chi2_muon/track.pid_n_dof);
    chi2_kaon->Fill(track.chi2_kaon/track.pid_n_dof);
    chi2_pion->Fill(track.chi2_pion/track.pid_n_dof);
    chi2_proton->Fill(track.chi2_proton/track.pid_n_dof);

    chi2_proton_m_muon->Fill(track.chi2_proton - track.chi2_muon);

  }
	  
  range_p->Fill(track.range_momentum); 
  mcs_p->Fill(track.mcs_momentum);
  deposited_e_max->Fill(track.deposited_energy_max);
  
  length->Fill(track.length);

  reco_momentum->Fill(track.momentum);
  is_contained->Fill(track.is_contained);

  dQdx_length->Fill(track.mean_trucated_dQdx, track.length);

  if (std::min(abs(track.start.Y() - border_y->GetBinLowEdge(1)), 
               abs(track.start.Y() - (border_y->GetBinLowEdge(border_y->GetNbinsX()) + border_y->GetBinWidth(border_y->GetNbinsX()))))
   <  std::min(abs(track.end.Y() - border_y->GetBinLowEdge(1)), 
               abs(track.end.Y() - (border_y->GetBinLowEdge(border_y->GetNbinsX()) + border_y->GetBinWidth(border_y->GetNbinsX()))))) {
    border_y->Fill(track.start.Y());
  }
  else {
    border_y->Fill(track.end.Y());
  }

  if (std::min(abs(track.start.X() - border_x->GetBinLowEdge(1)), 
               abs(track.start.X() - (border_x->GetBinLowEdge(border_x->GetNbinsX()) + border_x->GetBinWidth(border_x->GetNbinsX()))))
   <  std::min(abs(track.end.X() - border_x->GetBinLowEdge(1)), 
               abs(track.end.X() - (border_x->GetBinLowEdge(border_x->GetNbinsX()) + border_x->GetBinWidth(border_x->GetNbinsX()))))) {
    border_x->Fill(track.start.X());
  }
  else {
    border_x->Fill(track.end.X());
  }

  if (std::min(abs(track.start.Z() - border_z->GetBinLowEdge(1)), 
               abs(track.start.Z() - (border_z->GetBinLowEdge(border_z->GetNbinsX()) + border_z->GetBinWidth(border_z->GetNbinsX()))))
   <  std::min(abs(track.end.Z() - border_z->GetBinLowEdge(1)), 
               abs(track.end.Z() - (border_z->GetBinLowEdge(border_z->GetNbinsX()) + border_z->GetBinWidth(border_z->GetNbinsX()))))) {
    border_z->Fill(track.start.Z());
  }
  else {
    border_z->Fill(track.end.Z());
  }

  has_crt_hit_match->Fill(track.crt_match.hit.present);
  has_crt_track_match->Fill(track.crt_match.track.present);
  if (track.crt_match.track.present) {
    crt_match_time->Fill(track.crt_match.track.time);
    crt_track_angle->Fill(track.crt_match.track.angle);
  }
  else if (track.crt_match.hit.present) {
    crt_match_time->Fill(track.crt_match.hit.hit.time);
    crt_hit_distance->Fill(track.crt_match.hit.distance);
  }
  has_flash_match->Fill(track.flash_match.present);
  if (track.flash_match.present) {
    const numu::FlashMatch &flash_match = track.flash_match;
    double flash_time = flash_match.match_time_first;
    flash_match_time->Fill(flash_time);
    if (track.crt_match.track.present) {
      crt_v_flash_match_time->Fill(track.crt_match.track.time - flash_time);
    }
    else if (track.crt_match.hit.present) {
      crt_v_flash_match_time->Fill(track.crt_match.hit.hit.time - flash_time);
    }  
  }
  
  // check if truth match
  if (track.match.has_match && track.match.mcparticle_id >= 0) {
    const numu::RecoTrack &true_track = true_tracks.at(track.match.mcparticle_id);
    range_p_minus_truth->Fill(track.range_momentum - true_track.momentum);
    mcs_p_minus_truth->Fill(track.mcs_momentum - true_track.momentum); 
    deposited_e_max_minus_truth->Fill(track.deposited_energy_max - true_track.energy);
    deposited_e_avg_minus_truth->Fill(track.deposited_energy_avg - true_track.energy);
    deposited_e_med_minus_truth->Fill(track.deposited_energy_med - true_track.energy);
    
    range_p_diff->Fill(true_track.momentum, track.range_momentum - true_track.momentum);
    mcs_p_diff->Fill(true_track.momentum, track.mcs_momentum - true_track.momentum);
    deposited_e_max_diff->Fill(true_track.energy, track.deposited_energy_max - true_track.energy);
    
    range_p_comp->Fill(true_track.momentum, track.range_momentum);
    mcs_p_comp->Fill(true_track.momentum, track.mcs_momentum);
    deposited_e_max_comp->Fill(true_track.energy, track.deposited_energy_max);

    completion->Fill(track.match.completion);

    wall_enter->Fill(true_track.wall_enter);
    wall_exit->Fill(true_track.wall_exit);
    true_start_time->Fill(true_track.start_time);
  }
  else {
    completion->Fill(-0.5);
  }

  stopping_chisq_start->Fill(track.stopping_chisq_start);
  stopping_chisq_finish->Fill(track.stopping_chisq_finish);

  stopping_chisq->Fill(track.stopping_chisq_start);
  stopping_chisq->Fill(track.stopping_chisq_finish);

}
  } // namespace SBNOsc
} // namespace ana
