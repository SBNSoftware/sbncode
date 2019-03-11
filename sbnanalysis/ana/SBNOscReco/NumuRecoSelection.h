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
    mCC, mNC, mCosmic, mAll
  };

  std::string mode2Str(const InteractionMode &mode) const {
    switch (mode) {
      case mCC: return "CC";
      case mNC: return "NC";
      case mCosmic: return "Cosmic";
      case mAll: return "All";
    }
  }

  struct TrackInfo {
    double length;
    double energy;
    double open_angle;
    int pdgid;
    bool contained_in_cryo;
    bool contained_in_tpc;
    bool crosses_tpc;
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
    std::vector<std::vector<geo::BoxBoundedGeo>> tpc_volumes; //!< List of active tpc volumes -- retreived from Geoemtry service
    bool verbose; //!< Whether to print out info associated w/ selection.
    bool shakyMCTracks;
    std::vector<std::string> uniformWeights; //!< Weights taken from "EventWeight" that should be applied to the weight of each event
    double constantWeight; //!< constant weight to apply uniformly to each event
    double cosmicWeight;

    bool requireMatched;
    bool requireContained;

    std::string HitTag;
    std::string RecoTrackTag;
    std::string RecoVertexTag;
    std::string PFParticleTag;
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
  };

  static const unsigned nTruthCuts = 3; //!< number of truth cuts
  static const unsigned nRecoCuts = 3; //!< number of reco cuts
  static const unsigned nCuts = nTruthCuts + nRecoCuts;
  static const unsigned nModes = 4; //!< number of interaction modes
  static constexpr InteractionMode allModes[nModes] = {mCC, mNC, mCosmic, mAll};
  static constexpr const char* cutNames[nCuts] = {"Truth", "T_track", "T_match", "Reco", "R_match", "R_contained"};

  RecoEvent Reconstruct(const gallery::Event &ev, std::vector<RecoVertex> truth);
  std::vector<RecoParticle> RecoParticleInfo(const gallery::Event &event);
  std::vector<RecoParticle> SelectVertices(const std::vector<RecoParticle>& reco_particles);
  TrackInfo RecoTrackInfo(const std::vector<RecoParticle> &reco_particles, const RecoParticle& vertex);

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

  unsigned _event_counter;  //!< Count processed events
  unsigned _nu_count;  //!< Count selected events
  TGraph *_cut_counts; //!< Keep track of neutrinos per cut

  Config _config; //!< The config
  core::ProviderManager *_manager;

  RecoEvent _recoEvent;
  std::vector<RecoVertex> *_selected;

  RootHistos _root_histos[nCuts * nModes]; //!< Histos (one group per cut)

};

  }  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_NumuRecoSelection__

