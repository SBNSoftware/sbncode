#ifndef __sbnanalysis_HISTOGRAM_HH
#define __sbnanalysis_HISTOGRAM_HH

#include <string>
#include <vector>

#include <core/Event.hh>
#include "Cuts.h"

#include "../Data/RecoEvent.h"
#include "../Data/RecoTrack.h"
#include "../Data/Mode.h"

#include "DynamicSelector.h"

class TH1D;
class TH2D;
namespace ana {
 namespace SBNOsc {

/**
 *  Histograms associated with neutrino interactions. Filled for the list of all true
 *  and reco vertices. These histograms are constructed per interaction mode per cut.
 */
struct InteractionHistos {

  TH1D *track_length; //!< Length of the reconstructed primary track
  TH1D *nuE; //!< Neutrino energy
  TH1D *track_p; //!< Primary track momentum
  TH1D *true_deposited_energy; //!< Deposited energy of true track
  TH1D *beam_center_distance; //!< Distance of the neutrino interaction to the beam center
  TH1D *Q2; //!< Q2 of the interaction
  TH1D *true_contained_length; //!< True contained length of primary track
  TH1D *true_track_multiplicity; //!< True particle multiplicity of the interaction
  TH1D *crosses_tpc; //!< Whether the primary track crosses a TPC boundary
  TH1D *dist_to_match; //!< Distance from this vertex to the closest matching vertex reco->truth and truth->reco
  TH1D *primary_track_completion; //!< Completion of the primary track
  TH1D *n_reco_vertices; //!< Number of reconstructed vertices in the event with this vertex
  std::vector<TH1 *> all_histos;

  /**
 *  Intialize the histograms
 *  \param prefix A prefix to be added to each histogram name
 *  \param mode The mode of interaction for these histograms
 *  \param index The cut index for these histograms.
 */
  void Initialize(const std::string &prefix, numu::InteractionMode mode, unsigned index);
 
  /**
 * Scale all histograms by a value
 *
 * \param scale The amount to scale all histograms
 */
  void Scale(double scale);
  /**
 * Fill the histograms with a single interaction
 * \param vertex_index The index of this vertex into the list of truth/reco interactions
 * \param is_truth Whether this interaction is true or reco
 * \param event The reco event object
 * \param core_truth The list of true interactions from sbncode core
 */
  void Fill(
    unsigned vertex_index,
    bool is_truth,
    const numu::RecoEvent &event,
    const std::vector<event::Interaction> &core_truth);
  /**
 * Write the histograms to disk
 */
  void Write();
  /**
 * Add another set of histograms to this class
 * \param other The histograms to add to this class
 */
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
      case numu::mCCNonPrimary: return "CC-Other";
      case numu::mNC: return "NC";
      case numu::mNCNonPrimary: return "NC-Other";
      case numu::mCosmic: return "Cosmic";
      case numu::mOther: return "Other";
      case numu::mAll: return "All";
    }
  }
  static const unsigned nHistos = Cuts::nCuts + Cuts::nTruthCuts;
  static const unsigned nModes = 7; //!< number of interaction modes
  static constexpr numu::InteractionMode allModes[nModes] = 
    {numu::mCC, numu::mCCNonPrimary, numu::mNC, numu::mNCNonPrimary, numu::mCosmic, numu::mOther, numu::mAll}; //!< List of all interaction modes
  // static constexpr const char* histoNames[nHistos] = {"Truth", "Reco", "R_track", "R_vmatch", "R_tmatch", "R_match", "R_contained"}; //!< List of all cut names 
  static constexpr const char* histoNames[nHistos] = 
  {"Truth", "T_fid", "T_vqual", "T_tqual", "T_reco", 
   "Reco", "R_fid", "R_crttrack", "R_crthit", "R_length", "R_contained", "R_single", "R_crtactive"}; //!< Names of histograms
};

/**
 * Histograms to be filled per track
 */
struct TrackHistos {
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
  std::vector<TH1 *> all_histos;

  /**
 * Initialize this set of histograms
 * \param postfix The postfix to add to all histogram names
 */
  void Initialize(const std::string &postfix);
  /**
 * Scale all of the histograms in this class by a value
 * \param scale The scaling amount
 */
  void Scale(double scale);
  /**
 * Fill all of the histograms in this class with a track
 * \param track The track to fill
 * \param true_tracks The list of true particles in this event
 * \param cuts The configured cut class 
 */
  void Fill(
    const numu::RecoTrack &track,
    const std::map<size_t, numu::RecoTrack> &true_tracks,
    const Cuts &cuts);
  /**
 * Write this set of histograms to disk
 */
  void Write();
  /**
 * Add another set of histograms to this one
 * \param other The set of histograms to add
 */
  void Add(const TrackHistos &other);
  ~TrackHistos();

  /**
 * Check if a track is a certain mode for histograms
 * \param true_tracks Map of mcparticle-id's to true-track objects
 * \param track The reconstructed track to fill histograms with
 * \param mode_index The index of the mode in the range [0, nTrackHistos)
 * \return Whether the track in question is of the mode in question.
 */
  static bool IsMode(const std::map<size_t, numu::RecoTrack> &true_tracks, const numu::RecoTrack &track, unsigned mode_index);

  /**
 * Return the index of the track PDG into the list of PDG modes above
 * \param track The reconstructed track in question
 * \return The index in question 
 */
  static unsigned PDGIndex(const numu::RecoTrack &track);
};

/**
 * Set of Track and Interaction histos for all modes/cuts/sub-types
 *
 */
struct Histograms {

  InteractionHistos fInteraction[InteractionHistos::nHistos][InteractionHistos::nModes]; //!< all the interaction histograms
  std::vector<TrackHistos> fAllTracks; //!< Track histograms for all tracks
  std::vector<std::array<TrackHistos, Cuts::nCuts>> fPrimaryTracks; //!< Track histograms for priamry tracks in a candidate neutrino interaction

  void Initialize(const std::string &prefix="", std::vector<std::string> track_histo_types={});

  /**
 * Fill all of the histograms with an event
 * \param event The reconstructed event information
 * \param core The sbncode core event information
 * \param cuts The configured Cuts classes
 * \param fill_all_tracks Whether to fill all track histograms or just the primary track histograms
 */
  void Fill(const numu::RecoEvent &event, const event::Event &core, const Cuts &cuts, const std::vector<numu::TrackSelector> &selectors, bool fill_all_tracks=true);
  /**
 * Scale all histograms by a set value
 * \param scale The scaing value
 */
  void Scale(double scale);

  /**
 * Write all histograms to disk
 */
  void Write();
  /**
 * Add another set of histograms to this one
 * \param other The set of histograms to add
 */
  void Add(const Histograms &other);
};
 


  } // namespace SBNOSc
} // namespace ana

#endif
