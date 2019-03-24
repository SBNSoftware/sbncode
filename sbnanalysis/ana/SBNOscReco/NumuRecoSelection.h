#ifndef __sbnanalysis_ana_SBNOsc_NumuRecoSelection__
#define __sbnanalysis_ana_SBNOsc_NumuRecoSelection__

/**
 * \file NumuRecoSelection.h
 *
 * SBN nue selection.
 *
 * Author:
 */

#include <iostream>
#include <array>

#include "canvas/Utilities/InputTag.h"
#include "core/SelectionBase.hh"
#include "core/Event.hh"
#include "core/ProviderManager.hh"

#include "TH1D.h"
#include "TDatabasePDG.h"
#include "TGraph.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "larcorealg/Geometry/BoxBoundedGeo.h"

class TH2D;

namespace ana {
  namespace SBNOsc {

/**
 * \class NumuRecoSelection
 * \brief Electron neutrino event selection
 */
class NumuRecoSelection : public core::SelectionBase {
public:
  /** Constructor. */
  NumuRecoSelection();

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

  /**
 * Enum to hold each different typoe of reconstructed event
 */
  enum InteractionMode {
    mCC = 0, 
    mNC = 1, 
    mCosmic = 2, 
    mOther = 3,
    mAll = 4
  };

  /**
 *  Turn the InteractionMode enum into a string for (e.g.) histogram names.
 *
 *  \param mode The interaction mode to be converted
 *  \return String representaiton of that mode
 */
  std::string mode2Str(const InteractionMode &mode) const {
    switch (mode) {
      case mCC: return "CC";
      case mNC: return "NC";
      case mCosmic: return "Cosmic";
      case mOther: return "Other";
      case mAll: return "All";
    }
  }

  /**
 * Reconstruction information for the primary track of each muon neutrino vertex.
 * Produced for both reconstruction and truth information
 */
  struct TrackInfo {
    double length; //!< length of track
    double energy; //!< visible energy of track
    double costh; //!< cosine of angle to z axis
    int pdgid; //!< particle id
    bool contained_in_cryo; //!< is it contained a single cryostat?
    bool contained_in_tpc; //!< is it contained in a single TPC?
    bool crosses_tpc; //!< does it cross a tpc?
    bool is_contained; //!< is it contained in the "containment volume"?
    TVector3 start; //!< start position of track
    TVector3 end; //!< end position of track
  };

  /**
 *  Reconstruction information for each neutrino vertex.
 *  Produced from both reconstruction and truth information
 */
  struct RecoVertex {
    TVector3 position; //!< location of the vertex
    double nu_energy; //!< true/reconstructed neutrino energy
    TrackInfo track; //!< information on the primary track
    InteractionMode mode; //!< mode of the interaction
    int match; //!< index of the truth interaction matched to the reconstructed vertex. -1 if no match.
  };

  /** Reconstruction Information about Event */
  struct RecoEvent {
    std::vector<RecoVertex> reco; //!< List of reconstructed vertices
    std::vector<RecoVertex> truth; //!< List of truth vertices
  };

protected:
  /** Configuration parameters */
  struct Config {
    std::vector<geo::BoxBoundedGeo> fiducial_volumes; //!< List of FV containers -- set by "fiducial_volumes"
    std::vector<geo::BoxBoundedGeo> containment_volumes; //!< List of volumes for containment cuts -- set by "containment_volumes"
    std::vector<geo::BoxBoundedGeo> cryostat_volumes; //!< List of cryostat volumes -- retreived from Geometry service
    std::vector<std::vector<geo::BoxBoundedGeo>> tpc_volumes; //!< List of active tpc volumes -- retreived from Geoemtry service
    bool verbose; //!< Whether to print out info associated w/ selection.
    bool shakyMCTracks; //!< How to handle MC tracks with some missing truth information
    std::vector<std::string> uniformWeights; //!< Weights taken from "EventWeight" that should be applied to the weight of each event
    double constantWeight; //!< Constant weight to apply uniformly to each event
    double cosmicWeight; //!< Weight applied to all events matched to a cosmic track

    bool requireMatched; //!< Apply cut that requires each reconstructed vertex to be matched to a truth vertex
    bool requireTrack; //!< Apply cut that requires each reconstructed vertex to have an associated primary track
    bool requireContained; //!< Apply cut that requires each primary track to be contained inside the containment volume

    double cryostat_volume_boundary; 

    std::string HitTag; //!< art tag for hits
    std::string RecoTrackTag; //!< art tag for reconstructed tracks
    std::string RecoVertexTag; //!< art tag for reconstructed vertices
    std::string PFParticleTag; //!< art tag for PFParticles
    std::string CorsikaTag; //!< art tag for corsika MCTruth
    std::string CRTTrackTag; //!< art tag for CRT tracks
  };

  /** Histograms made for output 
 * Each histogram is made for each cut for each interaction mode
 *
 * */
  struct RootHistos {
    TH1D *track_length; //!< Length of the reconstructed primary track
  };

  /** Reconstructed information about each particle. Internal struct used
 * to combine information on each reconstructed particle. 
 * */
  struct RecoParticle {
    bool p_is_clear_cosmic; //!< Taken from Pandora metadata "is_clear_cosmic"
    bool p_is_neutrino; //!< Taken from Pandora metadata "is_neutrino"
    double p_nu_score; //!< Take from Pandora metadata "nu_score"
    std::vector<const recob::Vertex*> vertices; //!< List of vertices associated with the particle
    std::vector<size_t> daughters; //!< Daughters of the particle in the "particle flow". Value represents index into pandora information.
  };

  static const unsigned nTruthCuts = 3; //!< number of truth cuts
  static const unsigned nRecoCuts = 4; //!< number of reco cuts
  static const unsigned nCuts = nTruthCuts + nRecoCuts; //!< total number of cuts
  static const unsigned nModes = 5; //!< number of interaction modes
  static constexpr InteractionMode allModes[nModes] = {mCC, mNC, mCosmic, mOther, mAll}; //!< List of all interaction modes
  static constexpr const char* cutNames[nCuts] = {"Truth", "T_track", "T_match", "Reco", "R_track", "R_match", "R_contained"}; //!< List of all cut names 

  // Internal functions

  /**
 * Produce all reconstruction information from the gallery event
 * \param ev the gallery Event
 * \param truth the list of truth vertices for this event
 * \return the RecoEvent object for this event
 */
  RecoEvent Reconstruct(const gallery::Event &ev, std::vector<RecoVertex> truth);

  /**
 * Gathers together reconstruction information on each individual particle.
 * \param ev the gallery Event
 *
 * \return the list of RecoParticle objects for every reconstructed particle in the event
 */
  std::vector<RecoParticle> RecoParticleInfo(const gallery::Event &event);

  /**
 * Selects which RecoParticles are neutrino interaction candidates
 * \param reco_particles the list of all reconstructed particles to be considered
 *
 * \return the list of reconstructed particles which might be neutrinos
 */
  std::vector<RecoParticle> SelectVertices(const std::vector<RecoParticle>& reco_particles);

  /** 
 * Builds the reconstructed track information from a neutrino vertex candidate
 * \param ev the gallery Event
 * \param vertex the neutrino vertex candidate
 *
 * \return Track Information on the associated primary track. Returns nonsense if none is found
 */
  TrackInfo RecoTrackInfo(const gallery::Event &ev, const RecoParticle& vertex);

  /**
 *  Converts the NumuRecoSelection::RecoVertex information into reco information used by sbncode core
 *  \param truth the list of truth vertices produced by sbncode core
 *  \param vertex the reconstructed vertex produced by this selection
 *  \param weight the calculated weight for this interaction
 *
 *  \return Reconstruction information as needed by the sbncode core class
 */
  Event::RecoInteraction CoreRecoInteraction(const std::vector<Event::Interaction> &truth, const RecoVertex &vertex, double weight);

  /**
 * Produced vertex information from truth information
 * \param ev the gallery Event
 *
 * \return the list of truth neutrino interactions for this event 
 */
  std::vector<RecoVertex> MCTruthRecoVertexInfo(const gallery::Event &ev);

  /**
 * Get the primary track associated with a truth neutrino interaction.
 * \param ev the gallery Event
 * \param mctruth the MCTruth object for the considered truth neutrino interaction
 *
 * \return the index into the list of MCTrack's containing this track. Equal to -1 if no such track exists.
 */
  int GetNeutrinoTrack(const gallery::Event &ev, const simb::MCTruth &mctruth);

  /**
 * Get the primary track information associated with an mctruth object
 * \param ev the gallery Event
 * \param mc_truth the MCTruth object for the considered truth neutrino interaction
 *
 * \return information associated with the primary track. Set to nonsense if no such track exists.
 */
  TrackInfo MCTruthTrackInfo(const gallery::Event &event, const simb::MCTruth &mc_truth);

  /**
 * Get the TrackInfo information associated with an MCTrack
 * \param truth The neutrino interaction which produced this track
 * \param track The MCTrack information for this track
 *
 * \return reconstruction information associated with the truth track
 */
  TrackInfo MCTrackInfo(const simb::MCTruth &truth, const sim::MCTrack &track);

  /** 
 * Process each cut associated with reconstructed events
 * \param event The reconstructed event information
 * \param reco_vertex_index The index of the candidate reconstructed neutrino vertex into the list of such vertices in the RecoEvent
 *
 * \return A list of bool's of whether the reco event passes each cut
 */
  std::array<bool, nRecoCuts> ProcessRecoCuts(const RecoEvent &event, unsigned reco_vertex_index);

  /**
 * Select a reco event based on the cut values provided by ProcessRecoCuts
 * \param cuts the list of cuts returned by ProcessRecoCuts
 *
 * \return whether to select this reconstructed neutrino vertex candidate
 */
  bool SelectReco(std::array<bool, nRecoCuts> &cuts);

  /**
 * Process the cuts associated with each truth interaction
 * \param event The reconstructed event information
 * \param turth_vertex_index The index of the candidate truth neutrino vertex into the list of such vertices in the RecoEvent
 *
 * \return A list of bool's of whether the truth event passes each cut
 */
  std::array<bool, nTruthCuts> ProcessTruthCuts(const RecoEvent &event, unsigned truth_vertex_index);

  /** Helper function -- whether point is contained in fiducial volume list
 * \param v The point vector.
 *
 * \return Whether the point is contained in the configured list of fiducial volumes.
 * */
  bool containedInFV(const TVector3 &v);
  bool containedInFV(const geo::Point_t &v);

  /**
 * Calculate the associated topology of each primary track (i.e. whether the track is contained in things)
 * 
 * \param The reconstructed track information
 *
 * \return A list of bools associated with track topology. See the code for what is what.
 */
  std::array<bool, 4> RecoTrackTopology(const recob::Track *track);

  /**
 * Match a reconstructed vertex candidate to a cosmic track
 *
 * \param ev the gallery Event
 * \param vertex the reconstructed neutrino vertex candidate
 *
 * \return Whether that candidate can be matched to a cosmic track
 */
  bool MatchVertex2Cosmic(const gallery::Event &ev, const RecoVertex &vertex);

  unsigned _event_counter;  //!< Count processed events
  unsigned _nu_count;  //!< Count selected events
  TGraph *_cut_counts; //!< Keep track of neutrinos per cut

  Config _config; //!< The config
  core::ProviderManager *_manager; //!< Manager of larsoft services provided by sbncode core

  RecoEvent _recoEvent; //!< Branch container for the RecoEvent
  std::vector<RecoVertex> *_selected; //!< Branch container for the list of selected reco vertices

  RootHistos _root_histos[nCuts][nModes]; //!< Histos (one group per cut)

};

  }  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_NumuRecoSelection__

