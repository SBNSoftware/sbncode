#ifndef __sbnanalysis_ana_SBNOsc_NumuReco__
#define __sbnanalysis_ana_SBNOsc_NumuReco__

/**
 * \file NumuReco.h
 *
 * SBN nue selection.
 *
 * Author:
 */

#include <iostream>
#include <array>
#include <vector>

#include "canvas/Utilities/InputTag.h"
#include "core/SelectionBase.hh"
#include "core/Event.hh"
#include "core/ProviderManager.hh"

#include "TDatabasePDG.h"
#include "TGraph.h"

#include "canvas/Persistency/Common/FindManyP.h"
#include "lardataobj/Simulation/GeneratedParticleInfo.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

#include "larcorealg/Geometry/BoxBoundedGeo.h"

#include "LArReco/TrajectoryMCSFitter.h"
#include "LArReco/TrackMomentumCalculator.h"

#include "sbndcode/CRT/CRTProducts/CRTHit.hh"
#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"

#include "sbndcode/CRT/CRTUtils/CRTT0MatchAlg.h"
#include "sbndcode/CRT/CRTUtils/CRTTrackMatchAlg.h"

#include "../CosmicIDAlgs/ApaCrossCosmicIdAlg.h"
#include "../CosmicIDAlgs/StoppingParticleCosmicIdAlg.h"
#include "OpHitFinder/opHitFinderSBND.hh"


#include "../Data/CRTMatch.h"
#include "../Data/FlashMatch.h"
#include "../Data/Mode.h"
#include "../Data/RecoEvent.h"
#include "../Data/RecoParticle.h"
#include "../Data/RecoTrack.h"
#include "../Data/TruthMatch.h"
#include "../Data/MCType.h"

#include "../Histograms/CRTHisto.h"

class TH2D;
class TH1D;

namespace ana {
  namespace SBNOsc {

/**
 * \class NumuReco
 * \brief Electron neutrino event selection
 */
class NumuReco : public core::SelectionBase {
public:
  /** Constructor. */
  NumuReco();

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
  bool ProcessEvent(const gallery::Event& ev, const std::vector<event::Interaction> &truth, std::vector<event::RecoInteraction>& reco);



protected:
  /** Configuration parameters */
  struct Config {
    std::vector<geo::BoxBoundedGeo> containment_volumes; //!< List of volumes for containment cuts -- set by "containment_volumes"
    std::vector<geo::BoxBoundedGeo> active_volumes; //!< List of active volumes per cryostat
    std::vector<std::vector<geo::BoxBoundedGeo>> tpc_volumes; //!< List of active tpc volumes -- retreived from Geoemtry service

    bool verbose; //!< Whether to print out info associated w/ selection.
    std::vector<std::string> uniformWeights; //!< Weights taken from "EventWeight" that should be applied to the weight of each event
    double constantWeight; //!< Constant weight to apply uniformly to each event
    double cosmicWeight; //!< Weight applied to all events matched to a cosmic track

    bool requireMatched; //!< Apply cut that requires each reconstructed vertex to be matched to a truth vertex
    bool requireTrack; //!< Apply cut that requires each reconstructed vertex to have an associated primary track
    bool requireContained; //!< Apply cut that requires each primary track to be contained inside the containment volume
    double trackMatchContainmentCut;

    double CRTHitTimeCorrection;
    bool CRTHitinOpHitRange;
    double CRT2OPTimeWidth;

    bool CosmicIDAllTracks;

    bool MakeOpHits;

    int FlashMatchMethod;
    int TSMode;
    double flashMatchTimeDifference;

    unsigned PMTTriggerThreshold;

    std::array<float, 2> BeamSpillWindow;
    double beamCenterX;
    double beamCenterY;

    std::string RecoSliceTag;
    std::string RecoTrackTag; //!< art tag for reconstructed tracks
    std::string RecoVertexTag; //!< art tag for reconstructed vertices
    std::vector<std::string> TPCRecoTagSuffixes;
    std::string CaloTag;
    std::string PIDTag;
    std::string PFParticleTag; //!< art tag for PFParticles
    std::string FlashMatchTag;
    std::string CorsikaTag; //!< art tag for corsika MCTruth
    std::string CRTTrackTag; //!< art tag for CRT tracks
    std::string CRTHitTag;
    std::string OpFlashTag;
    std::string MCParticleTag; //!< art tag for MCParticle 

  };


  // Internal functions

  /**
 * Produce all reconstruction information from the gallery event
 * \param ev the gallery Event
 * \param truth the list of truth interactions for this event
 * \return the RecoEvent object for this event
 */
  numu::RecoEvent Reconstruct(const gallery::Event &ev, const std::vector<event::Interaction> &truth);

  /**
 *  Collect information on true particles for this event
 */
  std::map<size_t, numu::TrueParticle> MCParticleInfos();

  /**
 *  Adapt an MCParticle to the information needed for this event selection module
 */
  numu::TrueParticle MCParticleInfo(const simb::MCParticle &particle);

  /**
 * Return the list of intime CRTHits in the reconstructed event.
 * \return In time CRT Hits in the numu reco format
 */
  std::vector<numu::CRTHit> InTimeCRTHits();
  void FillCRTHits();

  /**
 *
 */
  numu::CRTHit SBND2numuCRTHit(const sbnd::crt::CRTHit &hit);

  /**
 * Returns whether the povided time is inside the configured beam spill.
 * \param time Input time in us.
 * \return Whether the provided time is inside the configured beam spill.
 */
  bool InBeamSpill(float time);

  /**
 * Gathers together reconstruction information on each individual particle.
 *
 * \return the list of RecoParticle objects for every reconstructed particle in the event
 */
  std::vector<numu::RecoParticle> RecoParticleInfo();

  /**
 *
 *  Returns whether this slice is designated as a neutrino interaction by Pandora
 *
 * \param slice The neutrino slice gathered from pandora information
 * \return Boolean which is true if the slice is a neutrino interaction
 */
  bool SelectSlice(const numu::RecoSlice &slice);

  /**
 * Gathers a map of track ID's to RecoTrack objects. This ID may not be the same as
 * the pandora ID.
 *
 * \return Map of track ID's to track objects
 */
  std::map<size_t, numu::RecoTrack> RecoTrackInfo();

  /**
 * Gathers the list of reconstructed candidate neutrino interactions as a list of
 * TPC slices.
 *
 * \param reco_tracks The list of reconstructed tracks in the event
 * \param particles The list of reconstructed PFParticles in the event
 */
  std::vector<numu::RecoSlice> RecoSliceInfo(
    std::map<size_t, numu::RecoTrack> &reco_tracks,
    const std::vector<numu::RecoParticle> &particles);

  /**
 * Returns whether a primary track (candidate muon) candidate exists for a 
 * neutrino interaction candidate
 * 
 * \param tracks The list of reconstructed tracks in the event
 * \param slice The candidate neutrino interaction
 * \return Whether a primary track candidate exists
 */
  bool HasPrimaryTrack(const std::map<size_t, numu::RecoTrack> &tracks, const numu::RecoSlice &slice);


  /**
 * Return the list of tracks associated with a neutrino interaction candidate.
 * \param tracks The list of all reconstructed tracks in the event
 * \param particles The list of reconstructed particles associated with this
 *                  neutrino interaction canidate 
 *
 * \return The list of track IDs associated with this slice
 */
  std::vector<size_t> RecoSliceTracks(
    const std::map<size_t, numu::RecoTrack> &tracks,
    const std::map<size_t, numu::RecoParticle> &particles);

  /**
 *  Converts the NumuReco::RecoInteraction information into reco information used by sbncode core
 *  \param truth the list of truth vertices produced by sbncode core
 *  \param vertex the reconstructed vertex produced by this selection
 *  \param weight the calculated weight for this interaction
 *
 *  \return Reconstruction information as needed by the sbncode core class
 */
  event::RecoInteraction CoreRecoInteraction(const std::vector<event::Interaction> &truth, const numu::RecoInteraction &vertex, double weight);

  /**
 *  Gets the ID of a true particle is stored by MC photon information
 *  \param mcparticle_id The MCParticleID of the true particle as returned by G4
 *
 *  \return The ID of the particle as in the photon objects
 */
  int GetPhotonMotherID(int mcparticle_id);
  /**
 *  Returns whether a TPC track has a match to Optical information
 *  \param pandora_track The LArSoft track object
 *  \param track The track infromation from this module
 *
 *  \return Flash Match object
 */
  numu::FlashMatch FlashMatching(const recob::Track &pandora_track, const numu::RecoTrack &track); 

  /**
 *  Returns whether a TPC track has a match to CRT information
 *  \param track The track information from this module
 *  \param pandora_track The LArSoft track object
 *  \param track_hits The list of hits associated with this track
 *
 *  \return CRTMatch object
 */
  numu::CRTMatch CRTMatching(const numu::RecoTrack &track, const recob::Track &pandora_track, const std::vector<art::Ptr<recob::Hit>> &track_hits);

  void ApplyCosmicID(numu::RecoTrack &track);

  /**
 * Calculate some topology factoids about a reconstructed track
 * 
 * \param track The pointer to LArSoft track information
 *
 * \return A list of bools associated with track topology. See the code for what is what.
 */
  std::array<bool, 4> RecoTrackTopology(const art::Ptr<recob::Track> &track);

  /**
 * Get matching information from a reconstructed track to truth information
 * \param track_id the "ID" of the reconstructed track. __NOT__ the pandora ID.
 *
 * \return The matchinf information of this track to truth information
 */
  numu::TrackTruthMatch MatchTrack2Truth(size_t track_id);

  /**
 * Fill up the Optical information containers from the gallery Event
 * \param ev The gallery Event
 */
  void CollectPMTInformation(const gallery::Event &ev);
  /**
 * Fill up the CRT information containers from the gallery Event
 * \param ev The gallery Event
 */
  void CollectCRTInformation(const gallery::Event &ev);
  /**
 * Fill up the TPC information containers from the gallery Event
 * \param ev The gallery Event
 */
  void CollectTPCInformation(const gallery::Event &ev);
  /**
 * Fill up truth information containers from the gallery Event
 * \param ev The gallery Event
 */
  void CollectTruthInformation(const gallery::Event &ev);

  /**
 * Test whether a point is in the configured fiducial volume
 * \param v The point to test
 *
 * \return Whether the point is in the configured fiducial volume
 */
  bool InActive(const TVector3 &v) const;

  unsigned _event_counter;  //!< Count processed events
  unsigned _nu_count;  //!< Count selected events
  TGraph *_cut_counts; //!< Keep track of neutrinos per cut

  Config _config; //!< The config

  // calculators for Reco things
  trkf::TrackMomentumCalculator *_track_momentum_calculator; //!< Calculator for range-based track momentum
  trkf::TrajectoryMCSFitter *_mcs_fitter; //!< Calculator for MCS based momentum
  

  numu::RecoEvent _recoEvent; //!< Branch container for the RecoEvent
  std::vector<numu::RecoInteraction> *_selected; //!< Branch container for the list of selected reco vertices

  sbnd::CRTTrackMatchAlg *_crt_track_matchalg; //!< Algorithm for matching reco Tracks -> CRT Tracks
  sbnd::CRTT0MatchAlg *_crt_hit_matchalg; //!< Algorithm for matching reco Tracks -> CRT hits (T0's)
  ApaCrossCosmicIdAlg _apa_cross_cosmic_alg; //!< Algorithm for doing cosmic ID by looking for tracks crossing APA
  StoppingParticleCosmicIdAlg _stopping_cosmic_alg; //!< Algorithm for doing cosmic ID using a fit to the energy deposits
  opdet::opHitFinderSBND *_op_hit_maker; //!< Optical hit maker

  // holders for CRT information
  const std::vector<sbnd::crt::CRTTrack> *_crt_tracks;
  std::vector<sbnd::crt::CRTTrack> _crt_tracks_local;
  std::vector<sbnd::crt::CRTHit> _crt_hits_local; 
  const std::vector<sbnd::crt::CRTHit> *_crt_hits; 
  bool _has_crt_hits;
  bool _has_crt_tracks;

  // holders for PMT information
  std::vector<art::Ptr<recob::OpHit>> _op_hit_ptrs;
  std::vector<recob::OpHit> _op_hits_local;

  // holders for TPC information
  std::vector<art::Ptr<recob::Slice>> _tpc_slices;
  std::vector<art::Ptr<recob::Track>> _tpc_tracks;
  std::vector<art::Ptr<recob::PFParticle>> _tpc_particles;
  std::vector<std::vector<art::Ptr<recob::PFParticle>>> _tpc_slices_to_particles;
  std::vector<std::vector<unsigned>> _tpc_slices_to_particle_index;
  std::vector<art::Ptr<recob::PFParticle>> _tpc_tracks_to_particles;
  std::vector<unsigned> _tpc_tracks_to_particle_index;
  std::map<unsigned, unsigned> _tpc_particles_to_track_index;
  std::vector<std::vector<art::Ptr<anab::Calorimetry>>> _tpc_tracks_to_calo;
  std::vector<std::vector<art::Ptr<anab::ParticleID>>> _tpc_tracks_to_pid;
  std::vector<std::vector<art::Ptr<recob::Hit>>> _tpc_tracks_to_hits;
  std::vector<std::vector<art::Ptr<anab::T0>>> _tpc_particles_to_T0;
  std::vector<std::vector<art::Ptr<recob::Vertex>>> _tpc_particles_to_vertex;
  std::vector<std::vector<unsigned>> _tpc_particles_to_daughters;
  std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> _tpc_particles_to_metadata;
  std::vector<std::vector<art::Ptr<anab::T0>>> _tpc_particles_to_flashT0;

  // holders for truth information
  std::vector<art::Ptr<simb::MCParticle>> _true_particles;
  std::map<int, art::Ptr<simb::MCTruth>> _true_particles_to_truth;
  std::map<int, const sim::GeneratedParticleInfo *> _true_particles_to_generator_info;

  // type of Monte Carlo
  numu::MCType fType;

  CRTHistos _crt_histograms;

};

  }  // namespace SBNOsc
}  // namespace ana

#endif  // __sbnanalysis_ana_SBNOsc_NumuReco__

