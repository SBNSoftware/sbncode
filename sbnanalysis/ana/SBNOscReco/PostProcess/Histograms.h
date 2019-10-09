#ifndef __sbnanalysis_HISTOGRAM_HH
#define __sbnanalysis_HISTOGRAM_HH

#include <string>
#include <vector>

#include <core/Event.hh>
#include "Cuts.h"

#include "../Data/RecoEvent.h"
#include "../Data/RecoTrack.h"
#include "../Data/Mode.h"

class TH1D;
class TH2D;
namespace ana {
 namespace SBNOsc {

struct InteractionHistos {
  TH1D *track_length; //!< Length of the reconstructed primary track
  TH1D *nuE;
  TH1D *track_p;
  TH1D *beam_center_distance;
  TH1D *Q2;
  TH1D *true_contained_length;
  TH1D *true_track_multiplicity;
  TH1D *crosses_tpc;
  TH1D *dist_to_match;
  TH1D *primary_track_completion;
  std::vector<TH1 *> all_histos;


  void Initialize(const std::string &prefix, numu::InteractionMode mode, unsigned index);
  void Scale(double scale);
  void Fill(
    unsigned vertex_index,
    bool is_truth,
    const numu::RecoEvent &event,
    const std::vector<event::Interaction> &core_truth);
  void Write();
  void Add(const InteractionHistos &other);
  ~InteractionHistos();
  
  /**
 *  Turn the InteractionMode enum into a string for (e.g.) histogram names.
 *
 *  \param mode The interaction mode to be converted
 *  \return String representaiton of that mode
 */
  std::string mode2Str(const numu::InteractionMode &mode) const {
    switch (mode) {
      case numu::mCC: return "CC";
      case numu::mNC: return "NC";
      case numu::mCosmic: return "Cosmic";
      case numu::mOther: return "Other";
      case numu::mAll: return "All";
    }
  }
  static const unsigned nHistos = Cuts::nCuts + Cuts::nTruthCuts;
  static const unsigned nModes = 5; //!< number of interaction modes
  static constexpr numu::InteractionMode allModes[nModes] = 
    {numu::mCC, numu::mNC, numu::mCosmic, numu::mOther, numu::mAll}; //!< List of all interaction modes
  // static constexpr const char* histoNames[nHistos] = {"Truth", "Reco", "R_track", "R_vmatch", "R_tmatch", "R_match", "R_contained"}; //!< List of all cut names 
  static constexpr const char* histoNames[nHistos] = {"Truth", "T_fid", "T_vqual", "T_tqual", "Reco", "R_fid", "R_passcrt", "R_contained"};
};

struct TrackHistos {
  TH1D *chi2_proton_diff;
  TH1D *chi2_muon_diff;
  TH1D *chi2_pion_diff;
  TH1D *chi2_kaon_diff;

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

  TH1D *has_crt_track_match;
  TH1D *has_crt_hit_match;
  TH1D *has_flash_match;

  TH1D *crt_match_time;
  TH1D *flash_match_time;
  TH1D *crt_v_flash_match_time;

  TH1D *stopping_chisq_start;
  TH1D *stopping_chisq_finish;
  std::vector<TH1 *> all_histos;

  static const unsigned nTrackHistos = 4;
  // static constexpr const char* trackHistoNames[nTrackHistos] = {"All", "Primary", "Contained", "Exiting", "Cosmic", "Neutrino", "No-Match"};
  static constexpr const char* trackHistoNames[nTrackHistos] = {"All", "Cosmic", "Neutrino", "No-Match"};

  void Initialize(const std::string &prefix, unsigned index);
  void Scale(double scale);
  void Fill(
    const numu::RecoTrack &track,
    const std::map<size_t, numu::RecoTrack> &true_tracks,
    const Cuts &cuts);
  void Write();
  void Add(const TrackHistos &other);
  ~TrackHistos();
};

struct Histograms {

  InteractionHistos fInteraction[InteractionHistos::nHistos][InteractionHistos::nModes];
  TrackHistos fAllTracks[TrackHistos::nTrackHistos];
  TrackHistos fPrimaryTracks[TrackHistos::nTrackHistos];

  Histograms(const std::string &prefix="");

  void Fill(const numu::RecoEvent &event, const event::Event &core, const Cuts &cuts, bool fill_all_tracks=true);
  void Scale(double scale);

  void Write();
  void Add(const Histograms &other);
};
 


  } // namespace SBNOSc
} // namespace ana

#endif
