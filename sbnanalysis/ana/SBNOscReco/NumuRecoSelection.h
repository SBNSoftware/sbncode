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

  enum InteractionMode {
    mCC = 0, 
    mNC = 1, 
    mCosmic = 2, 
    mOther = 3,
    mAll = 4
  };

  std::string mode2Str(const InteractionMode &mode) const {
    switch (mode) {
      case mCC: return "CC";
      case mNC: return "NC";
      case mCosmic: return "Cosmic";
      case mOther: return "Other";
      case mAll: return "All";
    }
  }

  struct TrackInfo {
    double length; // length of track
    double energy; // visible energy of track
    double costh; // cosine of angle to z axis
    int pdgid; // particle id
    bool contained_in_cryo; // is it contained a single cryostat?
    bool contained_in_tpc; // is it contained in a single TPC?
    bool crosses_tpc; // does it cross a tpc?
    bool is_contained;
    TVector3 start; // start position
    TVector3 end; // end position
  };

  struct RecoVertex {
    TVector3 position;
    double nu_energy;
    TrackInfo track;
    InteractionMode mode;
    int match;
  };

  /** Reconstruction Information about Event */
  struct RecoEvent {
    std::vector<RecoVertex> reco;
    std::vector<RecoVertex> truth;
  };

protected:
  /** Configuration parameters */
  struct Config {
    std::vector<geo::BoxBoundedGeo> fiducial_volumes; //!< List of FV containers -- set by "fiducial_volumes"
    std::vector<geo::BoxBoundedGeo> cryostat_volumes; //!< List of cryostat volumes -- retreived from Geometry service
    std::vector<geo::BoxBoundedGeo> containment_volumes;
    std::vector<std::vector<geo::BoxBoundedGeo>> tpc_volumes; //!< List of active tpc volumes -- retreived from Geoemtry service
    bool verbose; //!< Whether to print out info associated w/ selection.
    bool shakyMCTracks;
    std::vector<std::string> uniformWeights; //!< Weights taken from "EventWeight" that should be applied to the weight of each event
    double constantWeight; //!< constant weight to apply uniformly to each event
    double cosmicWeight;

    bool requireMatched;
    bool requireTrack;
    bool requireContained;

    double cryostat_volume_boundary;

    std::string HitTag;
    std::string RecoTrackTag;
    std::string RecoVertexTag;
    std::string PFParticleTag;
    std::string CorsikaTag;
    std::string CRTTrackTag;
  };

  /** Histograms made for output */
  struct RootHistos {
    TH1D *track_length;
  };

  /** Reconstructed information about each particle */
  struct RecoParticle {
    bool p_is_clear_cosmic;
    bool p_is_neutrino; 
    double p_nu_score;
    std::vector<const recob::Vertex*> vertices;
    std::vector<size_t> daughters;
  };

  static const unsigned nTruthCuts = 3; //!< number of truth cuts
  static const unsigned nRecoCuts = 4; //!< number of reco cuts
  static const unsigned nCuts = nTruthCuts + nRecoCuts;
  static const unsigned nModes = 5; //!< number of interaction modes
  static constexpr InteractionMode allModes[nModes] = {mCC, mNC, mCosmic, mOther, mAll};
  static constexpr const char* cutNames[nCuts] = {"Truth", "T_track", "T_match", "Reco", "R_track", "R_match", "R_contained"};

  RecoEvent Reconstruct(const gallery::Event &ev, std::vector<RecoVertex> truth);
  std::vector<RecoParticle> RecoParticleInfo(const gallery::Event &event);
  std::vector<RecoParticle> SelectVertices(const std::vector<RecoParticle>& reco_particles);
  TrackInfo RecoTrackInfo(const gallery::Event &ev, const RecoParticle& vertex);

  Event::RecoInteraction CoreRecoInteraction(const std::vector<Event::Interaction> &truth, const RecoVertex &vertex, double weight);

  std::vector<RecoVertex> MCTruthRecoVertexInfo(const gallery::Event &ev);
  int GetNeutrinoTrack(const gallery::Event &ev, const simb::MCTruth &mctruth);
  TrackInfo MCTruthTrackInfo(const gallery::Event &event, const simb::MCTruth &mc_truth);
  TrackInfo MCTrackInfo(const simb::MCTruth &truth, const sim::MCTrack &track);

  std::array<bool, nRecoCuts> ProcessRecoCuts(const RecoEvent &event, unsigned reco_vertex_index);
  bool SelectReco(std::array<bool, nRecoCuts> &cuts);
  std::array<bool, nTruthCuts> ProcessTruthCuts(const RecoEvent &event, unsigned truth_vertex_index);

  /** Helper function -- whether point is contained in fiducial volume list
 * \param v The point vector.
 *
 * \return Whether the point is contained in the configured list of fiducial volumes.
 * */
  bool containedInFV(const TVector3 &v);
  bool containedInFV(const geo::Point_t &v);

  // more helper functions
  std::array<bool, 4> RecoTrackTopology(const recob::Track *track);
  bool MatchVertex2Cosmic(const gallery::Event &ev, const RecoVertex &vertex);

  unsigned _event_counter;  //!< Count processed events
  unsigned _nu_count;  //!< Count selected events
  TGraph *_cut_counts; //!< Keep track of neutrinos per cut

  Config _config; //!< The config
  core::ProviderManager *_manager;

  RecoEvent _recoEvent;
  std::vector<RecoVertex> *_selected;

  RootHistos _root_histos[nCuts][nModes]; //!< Histos (one group per cut)

};

  }  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_NumuRecoSelection__

