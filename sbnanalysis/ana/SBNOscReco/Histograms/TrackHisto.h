#ifndef _sbnanalysis_TrackHisto_hh
#define _sbnanalysis_TrackHisto_hh

#include <vector>
#include <map>

#include "larcorealg/Geometry/BoxBoundedGeo.h"

#include "../Data/RecoEvent.h"
#include "HistoList.h"

class TH1D;
class TH2D;
class TProfile2D;
class TProfile;

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

  TH1D *n_daughters;

  TH1D *range_p;
  TH1D *mcs_p;

  TH1D *range_p_minus_truth;
  TH1D *mcs_p_minus_truth;

  TH2D *lengh_munus_truth_length;
  TH2D *range_p_minus_truth_length;
  TH2D *mcs_p_minus_truth_length;


  TH1D *length;
  TH1D *artlength;
  TH1D *lengthdiff;
  TH1D *reco_momentum;
  TH1D *is_contained;

  TH1D *completion; 
  TH1D *purity;

  TProfile *completion_x; 
  TProfile *purity_x;
 
  TH2D *range_p_diff;
  TH2D *mcs_p_diff;

  TH2D *range_p_comp;
  TH2D *mcs_p_comp;

  TH2D *dQdx_length;
  TH1D *border_y;
  TH1D *border_z;
  TH1D *border_x;
  TH1D *true_start_time;
  TH1D *true_start_time_zoom;

  TH1D *wall_enter;
  TH1D *wall_exit;

  TH1D *has_crt_track_match;
  TH1D *has_crt_hit_match;
  TH1D *crt_hit_distance;
  TH1D *crt_track_angle;

  TH1D *crt_match_time;

  TH1D *stopping_chisq_start;
  TH1D *stopping_chisq_finish;
  TH1D *stopping_chisq;

  TH1D *flash_match_time;
  TH1D *crt_v_flash_match_time;

  TH2D *pid_confusion_tr;
  TH1D *theta;
  TProfile2D *angular_chi2_proton;
  TProfile2D *angular_chi2_muon;
  TProfile2D *angular_pm_score;
  TProfile *chi2_muon_p;
  TProfile *chi2_proton_p;
  TProfile *pm_score_p;

  TH2D *collE_minus_truth_length;
  TH2D *bestplaneE_minus_truth_length;

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
    const std::map<size_t, numu::TrueParticle> &true_particles);

};
  } // namespace SBNOSc
} // namespace ana

#endif

