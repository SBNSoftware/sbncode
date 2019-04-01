#ifndef __sbnanalysis_ana_SBNOsc_NumuSelection__
#define __sbnanalysis_ana_SBNOsc_NumuSelection__

/**
 * \file NumuSelection.h
 *
 * SBN nue selection.
 *
 * Author:
 */

#include <iostream>
#include <array>
#include <map>

#include "canvas/Utilities/InputTag.h"
#include "core/SelectionBase.hh"
#include "core/Event.hh"

#include "Utilities.h"

#include "TH1D.h"
#include "TDatabasePDG.h"
#include "TGraph.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCTrack.h"

// take the geobox stuff from uboonecode
#include "ubcore/LLBasicTool/GeoAlgo/GeoAABox.h"
#include "ubcore/LLBasicTool/GeoAlgo/GeoAlgo.h"
#include "TRandom.h"
#include "TRandomGen.h"

class TH2D;

namespace ana {
  namespace SBNOsc {

/**
 * \class NumuSelection
 * \brief Electron neutrino event selection
 */
class NumuSelection : public core::SelectionBase {
public:
  /** Constructor. */
  NumuSelection();

  /**
   * Initialization.
   *
   * \param config A configuration, as a FHiCL ParameterSet object
   */
  void Initialize(fhicl::ParameterSet* config=NULL);

  /** Finalize and write objects to the output file. */
  void Finalize();

  /**
   * Process one event.
   *
   * \param ev A single event, as a gallery::Event
   * \param Reconstructed interactions
   * \return True to keep event
   */
  bool ProcessEvent(const gallery::Event& ev, const std::vector<Event::Interaction> &truth, std::vector<Event::RecoInteraction>& reco);

  /** Additional information used by the selection per neutrino interaction */
  struct NuMuInteraction {
    bool t_is_contained; //!< whether the (maybe faked) lepton track is totally contained in the active volume
    double t_contained_length; //!< the length of the (maybe faked) lepton track contained in the active volume [cm]
    double t_length; //!< total length of (maybe faked) lepton track [cm]
    int t_pdgid; //!< PDGID of primary track (muon or pi+)
    double t_energy_true; //!< True energy of primary track [GeV]
    double t_energy_smeared; //!< Smeared energy of primary track [GeV]
    TVector3 t_momentum;

    // default constructor -- fills with bogus info
    /*
    NuMuInteraction():
      t_is_contained(false),
      t_contained_length(-1),
      t_length(-1),
      t_pdgid(-1)
      {}*/
  };

protected:
  /** Configuration parameters */
  struct Config {
    bool doFVCut; //!< Whether to apply fiducial volume cut
    std::vector<geoalgo::AABox> fiducial_volumes; //!< List of FV containers -- set by "fiducial_volumes"
    std::vector<geoalgo::AABox> active_volumes; //!< List of active volumes
    double vertexDistanceCut; //!< Value of max distance [cm] between truth and reconstructed vertex. Will not apply cut if value is negative.
    bool verbose; //!< Whether to print out info associated w/ selection.
    bool acceptShakyTracks; //!< Whether to calculate length of tracks that don't have all points generated
    bool trajPointLength; //!< Whether to use trajectory points to calculate length
    double minLengthContainedTrack; //!< Minimum length [cm] of contained tracks. Will not apply cut if value is negative.
    double minLengthExitingTrack; //!< Minimum length [cm] of exiting tracks.  Will not apply cut if value is negative.
    double trackVisibleEnergyThreshold; //!< Energy threshold for track to be acounted in visible energy calculation [GeV].
    double showerEnergyDistortion; //!< Energy distortion of showers for visible energy calculation (%).
    double trackEnergyDistortion; //!< Energy distortion of tracks for visible energy calculation (%).
    double leptonEnergyDistortionContained; //<! Energy distortion of lepton (primary track) for visible energy calculation (%).
    double leptonEnergyDistortionLeavingA; //!< Parameter to be used in the energy distortion of primary lepton for visible energy calculation. 
      // (%) = -A * Log(B * L)  where L is the lepton contained length
    double leptonEnergyDistortionLeavingB; //!< Parameter to be used in the energy distortion of primary lepton for visible energy calculation. 
      // (%) = -A * Log(B * L)  where L is the lepton contained length
    bool cutKMEC; //!< Whether to remove MEC events (useful for studying difference w.r.t. proposal)
    bool onlyKMEC; //!< Whether to remove all non-MEC events
    int selectMode;
    int selectCCNC;
    double selectionEfficiency; //!< Signal efficiency weight applied to signal (charged current) events
    std::vector<std::string> uniformWeights; //!< Weights taken from "EventWeight" that should be applied to the weight of each event
    double constantWeight; //!< constant weight to apply uniformly to each event
  };

  /** Histograms made for output */
  struct RootHistos {
    TH1D *h_numu_ccqe; //!< histogram w/ CCQE energy veriable [GeV]
    TH1D *h_numu_trueE; //!< histogram w/ truth energy variable [GeV]
    TH1D *h_numu_visibleE; //!< histogram w/ visible energy variable (total muon momentum + kinetic hadron energy) [GeV]
    TH1D *h_numu_true_v_visibleE; //!< histogram w/ difference of visible and truth energy [GeV] 
    TH1D *h_numu_t_is_contained; //!< histogram w/ whether associated track is contained in AV
    TH1D *h_numu_contained_L; //!< histogram w/ FV contained length of track in CC event [cm]
    TH1D *h_numu_t_length; //!< histogram w/ total length of associated track [cm]
    TH1D *h_numu_t_is_muon; //!< histogram of whether associated track is a muon
    //TH1D *h_numu_m_energy_true; //!< histogram of truth energy of muon (excludes background pions)
    //TH1D *h_numu_m_energy_smeared; //!< histogram of smeared energy of muon (excludes background pions)
    TH2D *h_numu_Vxy; //!< 2D x-y vertex histogram [cm]
    TH2D *h_numu_Vxz; //!< 2D x-z vertex histogram [cm]
    TH2D *h_numu_Vyz; //!< 2D y-z vertex histogram [cm]
    TH1D *h_numu_Vx_sig; 
    TH1D *h_numu_Vy_sig; 
    TH1D *h_numu_Vz_sig; 
    TH1D *h_numu_Vx_bkg; 
    TH1D *h_numu_Vy_bkg; 
    TH1D *h_numu_Vz_bkg; 
    TH1D *h_numu_t_is_muon_sig; //!< histogram of whether associated track is a muon
    TH1D *h_numu_t_is_muon_bkg; //!< histogram of whether associated track is a muon
  };

  // helper struct holding track info -- see NumuInteraction for variable details
  struct TrackInfo {
    bool t_is_contained; //!< whether the (maybe faked) lepton track is totally contained in the active volume
    double t_contained_length; //!< the length of the (maybe faked) lepton track contained in the active volume [cm]
    double t_length; //!< total length of (maybe faked) lepton track [cm]
  };

  static const unsigned nCuts = 5; //!< number of cuts

  /* Applies FV cut
  * \param v The neutrino interaction vertex
  * \return Whether to apply FV cut on neutrino
  *
  * */
  bool passFV(const TVector3 &v) { return !_config.doFVCut ||containedInFV(v); }

  /** Applies reco-truth vertex matching cut
 * \param truth_v Truth vertex vector
 * \param reco_v Reconstructed vertex vector
 * \return Whether to apply reco-truth vertex matching cut: true == passed cut
 * */
  bool passRecoVertex(const TVector3 &truth_v, const TVector3 &reco_v);

  /** Applies truth length cut
 *  \param length Distance travelled by lepton
 *  \param stop_in_tpc Whether the lepton stopped in the TPC volume
 *  \return Whether to apply length cut: true == passed cut
 * */
  bool passMinLength(double length, bool stop_in_tpc);

  /** Run Selection on a neutrino
 *  \param ev The gallery event.
 *  \param mctruth The mctruth object associated with the currently considered interaction.
 *  \param truth_ind The index into the vector of MCTruth objects in the gallery event from which the "mctruth" object is
 *  \param intInfo The interaction info object associated with this mctruth object.
 *
 *  \return A list containing whether the interaction passed each cut (ordered in the same way as the "cutNames()")
 * */
  std::array<bool, nCuts> Select(const gallery::Event& ev, const simb::MCTruth& mctruth, unsigned truth_ind, const NumuSelection::NuMuInteraction &intInfo);

  /** Get associated interaction information from monte carlo
 * \param ev The gallery event.
 * \param mctruth The mctruth object associated with the currently considered interaction.
 *
 * \return NuMuInteraction object containing information from the mctruth object
 * */
  NuMuInteraction interactionInfo(const gallery::Event& ev, const simb::MCTruth &mctruth, VisibleEnergyCalculator &calculator);

 /** Get the interaction info associated with a track. This works right now because
 * all interaction info comes from the primary track. Might have to re-think structure
 * going forward.
 *
 * \param track The primary track.
 * \return Interaction info filled in from this track.
 */
  TrackInfo trackInfo(const sim::MCTrack &track);

  /** Helper function -- whether point is contained in fiducial volume list
 * \param v The point vector.
 *
 * \return Whether the point is contained in the configured list of fiducial volumes.
 * */
  bool containedInFV(const TVector3 &v);

  /** Helper function -- whether point is contained in active volume list 
 * \param v The point vector.
 *
 * \return Whether the point is contained in the configured list of active volumes.
 * */
  bool containedInAV(const TVector3 &v);

  unsigned _event_counter;  //!< Count processed events
  unsigned _nu_count;  //!< Count selected events
  TGraph *_cut_counts; //!< Keep track of neutrinos per cut

  /** Names of cuts
 * \return List of names of cuts (for histogram names)
 * */
  static const std::array<std::string, nCuts> cutNames() {
    return {"MEC", "AV", "Track", "FV", "min_L"};
  }

  TRandomMT64 _rand; //!< random number generation
  Config _config; //!< The config

  std::vector<NuMuInteraction> *_interactionInfo; //!< Branch holder

  RootHistos _root_histos[nCuts]; //!< Histos (one group per cut)

  std::map<std::string, double> _eventCategories;
};

  }  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_NumuSelection__

