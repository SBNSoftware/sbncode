#ifndef _sbnanalysis_TrackHisto_hh
#define _sbnanalysis_TrackHisto_hh

#include <vector>
#include <map>

#include "../Data/RecoEvent.h"
#include "HistoList.h"

class TH1D;
class TH2D;

namespace ana {
 namespace SBNOsc {
/**
 * Histograms to be filled per track
 */
struct TrackHistos : public HistoList {
  TH1D *chi2_proton_diff;
  TH1D *chi2_muon_diff;
  TH1D *chi2_pion_diff;
  TH1D *chi2_kaon_diff;

  TH1D *chi2_proton;
  TH1D *chi2_muon;
  TH1D *chi2_pion;
  TH1D *chi2_kaon;

  TH1D *chi2_proton_m_muon;

  TH1D *range_p;
  TH1D *mcs_p;
  TH1D *deposited_e_max;
  TH1D *deposited_e_avg;

  TH1D *range_p_minus_truth;
  TH1D *mcs_p_minus_truth;
  TH1D *deposited_e_max_minus_truth;
  TH1D *deposited_e_avg_minus_truth;
  TH1D *deposited_e_med_minus_truth;

  TH1D *length;
  TH1D *is_contained;

  TH1D *completion; 
 
  TH2D *range_p_diff;
  TH2D *mcs_p_diff;
  TH2D *deposited_e_max_diff;

  TH2D *range_p_comp;
  TH2D *mcs_p_comp;
  TH2D *deposited_e_max_comp;

  TH2D *dQdx_length;
  TH1D *border_y;
  TH1D *border_z;
  TH1D *true_start_time;

  TH1D *wall_enter;
  TH1D *wall_exit;

  TH1D *has_crt_track_match;
  TH1D *has_crt_hit_match;
  TH1D *has_flash_match;

  TH1D *crt_match_time;
  TH1D *flash_match_time;
  TH1D *crt_v_flash_match_time;

  TH1D *stopping_chisq_start;
  TH1D *stopping_chisq_finish;
  TH1D *stopping_chisq;

  /**
 * Initialize this set of histograms
 * \param postfix The postfix to add to all histogram names
 */
  void Initialize(const std::string &postfix);

  /**
 * Fill all of the histograms in this class with a track
 * \param track The track to fill
 * \param true_tracks The list of true particles in this event
 */
  void Fill(
    const numu::RecoTrack &track,
    const std::map<size_t, numu::RecoTrack> &true_tracks);
  ~TrackHistos();

};
  } // namespace SBNOSc
} // namespace ana

#endif

