#ifndef __sbnanalysis_HISTOGRAM_HH
#define __sbnanalysis_HISTOGRAM_HH

#include <vector>
#include <string>
#include <tuple>

#include "../Histograms/HistoList.h"
#include "../Histograms/TrackHisto.h"
#include "../Histograms/Profile.h"
#include "../Histograms/InteractionHisto.h"
#include "../Histograms/DynamicSelector.h"

#include "Cuts.h"

namespace ana {
  namespace SBNOsc {

/**
 * Set of Track and Interaction histos for all modes/cuts/sub-types
 */
struct Histograms : public HistoList {

  void Initialize(
  const std::string &prefix, 
  const std::vector<std::string> &track_histo_types, 
  const std::vector<std::string> &track_profile_types,
  const std::vector<std::tuple<unsigned, float, float>> &track_profile_xranges);

  /**
 * Fill all of the histograms with an event
 * \param event The reconstructed event information
 * \param core The sbncode core event information
 * \param cuts The configured Cuts classes
 * \param fill_all_tracks Whether to fill all track histograms or just the primary track histograms
 */
  void Fill(
  const numu::RecoEvent &event, 
  const event::Event &core, 
  const Cuts &cutmaker, 
  const std::vector<numu::TrackSelector> &selectors, 
  const std::vector<numu::ROOTValue> &xvalues, 
  bool fill_all_tracks=true);

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
      case numu::mIntimeCosmic: return "IntimeCosmic";
      case numu::mOther: return "Other";
      case numu::mAll: return "All";
    }
  }
  static const unsigned nHistos = Cuts::nCuts + Cuts::nTruthCuts;
  static const unsigned nModes = 8; //!< number of interaction modes
  static constexpr numu::InteractionMode allModes[nModes] = 
    {numu::mCC, numu::mCCNonPrimary, numu::mNC, numu::mNCNonPrimary, numu::mCosmic, numu::mIntimeCosmic, numu::mOther, numu::mAll}; //!< List of all interaction modes
  // static constexpr const char* histoNames[nHistos] = {"Truth", "Reco", "R_track", "R_vmatch", "R_tmatch", "R_match", "R_contained"}; //!< List of all cut names 
  static constexpr const char* histoNames[nHistos] = 
  {"Truth", "T_fid", "T_vqual", "T_tqual", "T_reco", 
   "Reco", "R_fid", "R_goodmcs", "R_crttrack", "R_crthit", "R_length", "R_contained", "R_crtactive"}; //!< Names of histograms

  InteractionHistos fInteraction[nHistos][nModes]; //!< all the interaction histograms
  std::vector<TrackHistos> fAllTracks; //!< Track histograms for all tracks
  std::vector<std::array<TrackHistos, Cuts::nCuts>> fPrimaryTracks; //!< Track histograms for priamry tracks in a candidate neutrino interaction
  std::vector<std::vector<std::array<TrackProfiles, Cuts::nCuts>>> fPrimaryTrackProfiles; //!< Profile histograms for primary tracks

};

  } // namespace SBNOSc
} // namespace ana

#endif
