#ifndef _sbnanalysis_TrackHisto_hh
#define _sbnanalysis_TrackHisto_hh

#include <vector>
#include <map>

#include "larcorealg/Geometry/BoxBoundedGeo.h"

#include "../Data/RecoEvent.h"
#include "../MultiThread/THShared.h"
#include "HistoList.h"

class TH1D;
class TH2D;

namespace ana {
 namespace SBNOsc {
/**
 * Histograms to be filled per track
 */
struct TrackHistos : public HistoList {
  TH1Shared chi2_proton_diff;
  TH1Shared chi2_muon_diff;
  TH1Shared chi2_pion_diff;
  TH1Shared chi2_kaon_diff;

  TH1Shared chi2_proton;
  TH1Shared chi2_muon;
  TH1Shared chi2_pion;
  TH1Shared chi2_kaon;

  TH1Shared chi2_proton_m_muon;

  TH1Shared range_p;
  TH1Shared mcs_p;
  TH1Shared deposited_e_max;
  TH1Shared deposited_e_avg;

  TH1Shared range_p_minus_truth;
  TH1Shared mcs_p_minus_truth;

  TH2Shared range_p_minus_truth_length;
  TH2Shared mcs_p_minus_truth_length;

  TH1Shared deposited_e_max_minus_truth;
  TH1Shared deposited_e_avg_minus_truth;
  TH1Shared deposited_e_med_minus_truth;

  TH1Shared length;
  TH1Shared artlength;
  TH1Shared lengthdiff;
  TH1Shared reco_momentum;
  TH1Shared is_contained;

  TH1Shared completion; 
 
  TH2Shared range_p_diff;
  TH2Shared mcs_p_diff;
  TH2Shared deposited_e_max_diff;

  TH2Shared range_p_comp;
  TH2Shared mcs_p_comp;
  TH2Shared deposited_e_max_comp;

  TH2Shared dQdx_length;
  TH1Shared border_y;
  TH1Shared border_z;
  TH1Shared border_x;
  TH1Shared true_start_time;
  TH1Shared true_start_time_zoom;

  TH1Shared wall_enter;
  TH1Shared wall_exit;

  TH1Shared has_crt_track_match;
  TH1Shared has_crt_hit_match;
  TH1Shared crt_hit_distance;
  TH1Shared crt_track_angle;

  TH1Shared crt_match_time;

  TH1Shared stopping_chisq_start;
  TH1Shared stopping_chisq_finish;
  TH1Shared stopping_chisq;

  TH1Shared flash_match_time;
  TH1Shared crt_v_flash_match_time;

  TH2Shared pid_confusion_tr;

  /**
 * Initialize this set of histograms
 * \param postfix The postfix to add to all histogram names
 */
  void Initialize(const std::string &postfix, const geo::BoxBoundedGeo &detector_volume, double max_length);

  /**
 * Fill all of the histograms in this class with a track
 * \param track The track to fill
 * \param true_tracks The list of true particles in this event
 */
  void Fill(
    const numu::RecoTrack &track,
    const std::map<size_t, numu::RecoTrack> &true_tracks);

};
  } // namespace SBNOSc
} // namespace ana

#endif

