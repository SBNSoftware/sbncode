#ifndef __sbnanalysis_HISTOGRAM_HH
#define __sbnanalysis_HISTOGRAM_HH

#include <vector>
#include <string>
#include <tuple>

#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larcorealg/Geometry/GeometryCore.h"

#include "sbndcode/Geometry/GeometryWrappers/CRTGeoAlg.h"

#include "../Histograms/HistoList.h"
#include "../Histograms/TrackHisto.h"
#include "../Histograms/TrueParticleHisto.h"
#include "../Histograms/Profile.h"
#include "../Histograms/InteractionHisto.h"
#include "../Histograms/CosmicHisto.h"
#include "../Histograms/CRTHisto.h"
#include "../Histograms/DynamicSelector.h"

#include "Cuts.h"

namespace ana {
  namespace SBNOsc {

  // boilerplate to constexpr-append two arrays
  template <typename T, std::size_t N1, std::size_t N2>
  constexpr std::array<T, N1+N2> concat(std::array<T, N1> lhs, std::array<T, N2> rhs) {
    std::array<T,N1+N2> result{};
    size_t index = 0;
    for (auto &el : lhs)  {
      result[index] = std::move(el);
      index ++;
    }
    for (auto &el : rhs) {
      result[index] = std::move(el);
      index ++;
    }
    return result;
  }
  // std::array<const char *, Cuts::nTruthCuts+Cuts::nCuts> concat<const char *, Cuts::nTruthCuts, Cuts::nCuts>();

/**
 * Set of Track and Interaction histos for all modes/cuts/sub-types
 */
struct Histograms : public HistoList {

  void Initialize(
  const geo::GeometryCore *geometry,
  const sbnd::CRTGeoAlg &crt_geo,
  const Cuts &cuts,
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
  const std::vector<numu::TrackFunction> &xfunctions,
  bool fill_all_tracks=true);

  /**
 *  Turn the InteractionMode enum into a string for (e.g.) histogram names.
 *
 *  \param mode The interaction mode to be converted
 *  \return String representaiton of that mode
 */
  static std::string mode2Str(const numu::InteractionMode &mode) {
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
    { numu::mCC, numu::mCCNonPrimary, 
      numu::mNC, numu::mNCNonPrimary, 
      numu::mCosmic, numu::mIntimeCosmic, 
      numu::mOther, numu::mAll }; //!< List of all interaction modes

  static const unsigned nPIDs = 5;
  static constexpr const char * allPIDs[nPIDs] = 
  { "mu", "pi", "p", "k", "other"};

  // static constexpr const char* histoNames[nHistos] = {"Truth", "Reco", "R_track", "R_vmatch", "R_tmatch", "R_match", "R_contained"}; //!< List of all cut names 


  static constexpr std::array<const char *, nHistos> histoNames = concat(Cuts::truthCutNames, Cuts::cutNames); //!< Names of histograms

  InteractionHistos fInteraction[nHistos][nModes]; //!< all the interaction histograms
  std::vector<TrackHistos> fAllTracks; //!< Track histograms for all tracks
  std::vector<std::array<TrackHistos, Cuts::nCuts>> fPrimaryTracks; //!< Track histograms for priamry tracks in a candidate neutrino interaction
  std::vector<std::vector<std::array<TrackProfiles, Cuts::nCuts>>> fPrimaryTrackProfiles; //!< Profile histograms for primary tracks
  std::array<CosmicHistos, 4> fCosmic;
  TrueParticleHistos fParticles[nPIDs][2];

  std::array<CRTHistos, Cuts::nCuts> fCRTs;

};

  } // namespace SBNOSc
} // namespace ana

#endif
