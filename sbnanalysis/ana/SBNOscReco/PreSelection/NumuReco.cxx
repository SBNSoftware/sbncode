#include <list>
#include <string>
#include <iostream>
#include <cmath>
#include <vector>
#include <queue>

#include "fhiclcpp/ParameterSet.h"

#include <TH2D.h>
#include <TH1D.h>

#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "lardataobj/RecoBase/OpHit.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/ParticleID.h"

#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/GeometryCore.h"

#include "core/Event.hh"
#include "core/Experiment.hh"
#include "NumuReco.h"
#include "../SBNOsc/Utilities.h"
#include "../RecoUtils/RecoUtils.h"
#include "../RecoUtils/GeoUtil.h"

#include "ubcore/LLBasicTool/GeoAlgo/GeoAABox.h"
#include "canvas/Persistency/Provenance/ProductID.h"

#include "larsim/MCCheater/ParticleInventory.h"
#include "larsim/MCCheater/PhotonBackTracker.h"
#include "icaruscode/CRT/CRTProducts/CRTHit.hh"
#include "icaruscode/CRT/CRTProducts/CRTTrack.hh"
#include "sbndcode/CRT/CRTUtils/TPCGeoUtil.h"

// copied in RecoUtils here to not use art services

namespace ana {
  namespace SBNOsc {

const static TVector3 InvalidTVector3 = TVector3(-999, -999, -999);

void DumpTrueStart(const gallery::Event &ev, int mcparticle_id) {
  // track.match.mcparticle_id);
  const std::vector<simb::MCParticle> &mcparticle_list = *ev.getValidHandle<std::vector<simb::MCParticle>>("largeant");
      std::cout << "Track: " << mcparticle_id;
  for (const simb::MCParticle &part: mcparticle_list) {
    if (mcparticle_id == part.TrackId()) {
      std::cout << " start T: " << part.Position().T() << " to: " << part.EndPosition().T() << std::endl;
    }
  }
  return;
}

int NumuReco::GetPhotonMotherID(int mcparticle_id) {
  int index = -1;
  for (unsigned i = 0; i < _true_particles.size(); i++) {
    if (_true_particles[i]->TrackId() == mcparticle_id) {
      index = (int)i;
      break;
    }
  }
  assert(index >= 0);
  bool is_cosmic = _true_particles_to_truth.at(index)->Origin() != simb::Origin_t::kBeamNeutrino;
  if (is_cosmic) {
    unsigned truth_index = _true_particles_to_generator_info.at(index)->generatedParticleIndex();
    return (int) truth_index +1;
  }
  else {
    return mcparticle_id;
  }
}

sbnd::crt::CRTHit ICARUS2SBNDCrtHit(const icarus::crt::CRTHit& inp) {
  sbnd::crt::CRTHit ret;
  ret.feb_id = inp.feb_id;
  ret.pesmap = inp.pesmap;
  ret.peshit = inp.peshit;
  ret.ts0_s = inp.ts0_s;
  ret.ts0_s_corr = inp.ts0_s_corr;
  ret.ts0_ns = inp.ts0_ns;
  ret.ts0_ns_corr = inp.ts0_ns_corr;
  ret.ts1_ns = inp.ts1_ns;
  ret.plane = inp.plane;
  ret.x_pos = inp.x_pos;
  ret.x_err = inp.x_err;
  ret.y_pos = inp.y_pos;
  ret.y_err = inp.y_err;
  ret.z_pos = inp.z_pos;
  ret.z_err = inp.z_err;
  ret.tagger = inp.tagger;
  return ret;
}

NumuReco::NumuReco() :
  SelectionBase(),
  _event_counter(0),
  _nu_count(0),
  _selected(new std::vector<numu::RecoInteraction>) {}

void NumuReco::Initialize(fhicl::ParameterSet* config) {
  if (config) {
    fhicl::ParameterSet pconfig = config->get<fhicl::ParameterSet>("NumuReco");

    // configure the Cosmic ID algorithms
    _crt_track_matchalg = new sbnd::CRTTrackMatchAlg(fhicl::Table<sbnd::CRTTrackMatchAlg::Config>(pconfig.get<fhicl::ParameterSet>("CRTTrackMatchAlg"))(), 
      fProviderManager->GetGeometryProvider(), fProviderManager->GetDetectorPropertiesProvider());

    _crt_hit_matchalg = new sbnd::CRTT0MatchAlg(fhicl::Table<sbnd::CRTT0MatchAlg::Config>(pconfig.get<fhicl::ParameterSet>("CRTT0MatchAlg"))(),
      fProviderManager->GetGeometryProvider(), fProviderManager->GetDetectorPropertiesProvider());

    _apa_cross_cosmic_alg.reconfigure(*fProviderManager, fhicl::Table<ApaCrossCosmicIdAlg::Config>(pconfig.get<fhicl::ParameterSet>("ApaCrossCosmicIdAlg"))());

    _stopping_cosmic_alg.reconfigure(*fProviderManager, fhicl::Table<StoppingParticleCosmicIdAlg::Config>(pconfig.get<fhicl::ParameterSet>("StoppingParticleCosmicIdAlg"))());


    _config.tpc_volumes = SBNRecoUtils::TPCVolumes(fProviderManager->GetGeometryProvider());
    _config.active_volumes = SBNRecoUtils::ActiveVolumes(fProviderManager->GetGeometryProvider());

    // get the beam center
    _config.beamCenterX = pconfig.get<float>("beamCenterX", 130.);
    _config.beamCenterY = pconfig.get<float>("beamCenterY", 0.);

    _config.shakyMCTracks = pconfig.get<bool>("shakyMCTracks", false);
    _config.verbose = pconfig.get<bool>("verbose", false);

    _config.requireTrack = pconfig.get<bool>("requireTrack", false);

    _config.trackMatchContainmentCut = pconfig.get<double>("trackMatchContainmentCut", -1);

    _config.TSMode = pconfig.get<int>("TSMode", 1);
    _config.MakeOpHits = pconfig.get<bool>("MakeOpHits", false);
    _config.CRTHitDist = pconfig.get<double>("CRTHitDist", -1);
    
    _config.requireMatched = pconfig.get<bool>("requireMatched", false);
    _config.requireContained = pconfig.get<bool>("requireContained", false);

    // setup weight config
    _config.uniformWeights = pconfig.get<std::vector<std::string>>("uniformWeights", {});
    _config.constantWeight = pconfig.get<double>("constantWeight", 1.0);
    _config.cosmicWeight = pconfig.get<double>("cosmicWeight", 1.0);

    // flash match method
    _config.FlashMatchMethod = pconfig.get<int>("FlashMatchMethod", 2);
    _config.flashMatchTimeDifference = pconfig.get<double>("flashMatchTimeDifference");


    _config.CosmicIDAllTracks = pconfig.get<bool>("CosmicIDAllTracks", false);

    // whether to use flash matching in makin CRT decision
    _config.CRT2OPTimeWidth = pconfig.get<double>("CRT2OPTimeWidth", 0.);
    _config.CRTHitinOpHitRange = pconfig.get<bool>("CRTHitinOpHitRange", false);

    // get tag names
    _config.RecoTrackTag = config->get<std::string>("RecoTrackTag", "pandoraTrack");
    _config.RecoSliceTag = config->get<std::string>("RecoSliceTag", "pandora");
    _config.RecoVertexTag = config->get<std::string>("RecoVertexTag", "pandora");
    _config.PFParticleTag = config->get<std::string>("PFParticleTag", "pandora");
    _config.CaloTag = config->get<std::string>("CaloTag", "pandoraCalo");
    _config.PIDTag = config->get<std::string>("PIDTag", "pandoraPid");
    _config.CorsikaTag = config->get<std::string>("CorsikaTag", "cosmgen");
    _config.CRTTrackTag = config->get<std::string>("CRTTrackTag", "crttrack");
    _config.CRTHitTag = config->get<std::string>("CRTHitTag", "crthit");
    _config.OpFlashTag = config->get<std::string>("OpFlashTag", "ophit");
    _config.MCParticleTag = config->get<std::string>("MCParticleTag", "largeant");
    _config.TPCRecoTagSuffixes = config->get<std::vector<std::string>>("TPCRecoTagSuffixes", { "" });

    {
      fhicl::ParameterSet dCV = \
        pconfig.get<fhicl::ParameterSet>("containment_volume_inset");
      double dx = dCV.get<double>("x");
      double dy = dCV.get<double>("y");
      double zfront = dCV.get<double>("zfront");
      double zback = dCV.get<double>("zback");
      for (const geo::BoxBoundedGeo &geo: _config.active_volumes) {
        _config.containment_volumes.emplace_back(geo.MinX() + dx, geo.MaxX() - dx, geo.MinY() + dy, geo.MaxY() - dy, geo.MinZ() + zfront, geo.MaxZ() - zback);
      }
    }

    // initialize things to do with reconstruction
    double min_track_length = pconfig.get<double>("MinTrackLength", 10.);
    _track_momentum_calculator = new trkf::TrackMomentumCalculator(min_track_length);
    _mcs_fitter = new trkf::TrajectoryMCSFitter(pconfig.get<fhicl::ParameterSet>("MCSFitter", {}));

    _op_hit_maker = (_config.MakeOpHits) ? new opdet::opHitFinderSBND(pconfig.get<fhicl::ParameterSet>("OpHitMaker"), fProviderManager->GetDetectorClocksProvider()) :
        NULL;
  }

  // Setup histo's for root output
  fOutputFile->cd();

  // add branches
  fTree->Branch("reco_event", &_recoEvent);
  fTree->Branch("reco_vertices", &_selected);


  hello();
}


void NumuReco::Finalize() {
 // cleanup pointers to algorithms
 if (_crt_track_matchalg != NULL) delete _crt_track_matchalg;
 if (_crt_hit_matchalg != NULL) delete _crt_hit_matchalg;
 if (_op_hit_maker != NULL) delete _op_hit_maker;
 if (_track_momentum_calculator != NULL) delete _track_momentum_calculator;
 if (_mcs_fitter != NULL) delete _mcs_fitter;
}

double NumuReco::TrackCompletion(int mcparticle_id, const std::vector<art::Ptr<recob::Hit>> &reco_track_hits) {
  // get handle to back tracker
  cheat::BackTracker *bt = fProviderManager->GetBackTrackerProvider();

  // get all the IDE's of the truth track
  const std::vector<const sim::IDE*> mcparticle_ides = bt->TrackIdToSimIDEs_Ps(mcparticle_id);
  // sum it up
  double mcparticle_energy = 0.;
  for (auto const &ide: mcparticle_ides) {
    mcparticle_energy += ide->energy;
  }

  // get all the hits of the reco track that match the truth track
  const std::vector<art::Ptr<recob::Hit>> matched_reco_track_hits = bt->TrackIdToHits_Ps(mcparticle_id, reco_track_hits);

  // for each of the hits get the energy coming from the track
  double matched_reco_energy = 0.;
  for (auto const &matched_reco_track_hit: matched_reco_track_hits) {
    std::vector<sim::IDE> this_hit_IDEs = bt->HitToAvgSimIDEs(*matched_reco_track_hit);
    for (auto const &ide: this_hit_IDEs) {
      if (ide.trackID == mcparticle_id) {
        matched_reco_energy += ide.energy;
      }
    }
  }

  return matched_reco_energy / mcparticle_energy;
}

event::RecoInteraction NumuReco::CoreRecoInteraction(const std::vector<event::Interaction> &truth, const numu::RecoInteraction &vertex, double weight) {
  event::RecoInteraction ret;
  if (vertex.match.mctruth_vertex_id >= 0) {
    ret.truth_index = vertex.match.mctruth_vertex_id;
  }
  ret.reco_energy = vertex.nu_energy;
  ret.weight = weight;
  return ret;
}

bool NumuReco::ProcessEvent(const gallery::Event& ev, const std::vector<event::Interaction> &core_truth, std::vector<event::RecoInteraction>& reco) {
  if (_event_counter % 10 == 0) {
    std::cout << "NumuReco: Processing event " << _event_counter << " "
              << "(" << _nu_count << " neutrinos selected)"
              << std::endl;
  }
  std::cout << "New Event! " << _event_counter << std::endl;
  // clear out old containers
  _selected->clear();

  bool selected = false;

  _event_counter++;

  std::map<size_t, numu::RecoTrack> true_tracks = MCParticleTracks(ev);

  // get the truth vertex info
  std::vector<numu::RecoInteraction> truth = MCTruthInteractions(ev, true_tracks);

  // get the event info
  _recoEvent = Reconstruct(ev, truth);

  // complete the information
  _recoEvent.true_tracks = std::move(true_tracks);

  // save the information
  for (unsigned i = 0; i < _recoEvent.reco.size(); i++) {
    const numu::RecoInteraction &vertex = _recoEvent.reco[i];
    // compute the weight for this interaction
    double weight = 1.;
    weight *= _config.constantWeight;
    // TODO: what about cosmics?
    // if (vertex.match.mctruth_vertex_id >= 0) {
    //   for (auto const &key: _config.uniformWeights) {
    //     weight *= core_truth[vertex.match.mctruth_vertex_id].weightmap.at(key)[0];
    //  }
    // }
    if (vertex.match.mode == numu::mCosmic) {
      weight *= _config.cosmicWeight;
    }

    // store reco interaction
    _selected->push_back(vertex);
    // store sbncode reco interaction
    reco.push_back(CoreRecoInteraction(core_truth, vertex, weight));
    selected = true;
    _nu_count++;
  }

  return true;
}

numu::TrackTruthMatch MCTrackMatch(const simb::MCTruth &truth, const simb::MCParticle &track) {
  // get the match info
  numu::TrackTruthMatch track_match;
  track_match.has_match = true;
  track_match.mctruth_has_neutrino = truth.NeutrinoSet();
  track_match.mctruth_vertex = truth.NeutrinoSet() ? truth.GetNeutrino().Nu().Position().Vect() : InvalidTVector3;
  track_match.mctruth_origin = truth.Origin();
  track_match.mctruth_ccnc = truth.NeutrinoSet() ? truth.GetNeutrino().CCNC() : -1;
  track_match.mcparticle_id = track.TrackId();
  track_match.completion = 1.;
  track_match.match_pdg = track.PdgCode();
  return track_match;
}

// get information associated with track
numu::RecoTrack NumuReco::MCTrackInfo(const simb::MCParticle &track) {
  // default values
  double contained_length = 0;
  bool crosses_tpc = false;

  // If the interaction is outside the active volume, then g4 won't generate positions for the track.
  // So size == 0 => outside FV
  //
  // If size != 0, then we have to check volumes
  bool contained_in_cryo = track.NumberTrajectoryPoints() > 0;
  bool contained_in_tpc = track.NumberTrajectoryPoints() > 0;


  bool is_contained = true;

  // other truth information
  double costh = track.Pz() / track.Momentum().Vect().Mag();
  int pdgid = track.PdgCode();
  double kinetic_energy = track.E() /* already in GeV*/ - PDGMass(pdgid) / 1000. /* MeV -> GeV */;

  // setup intial track locations
  TLorentzVector pos = track.Position();
  // get the active volume that the start position is in
  int cryostat_index = -1;
  int tpc_index = -1;
  int containment_index = -1;
  for (int i = 0; i < _config.containment_volumes.size(); i++) {
    if (_config.containment_volumes[i].ContainsPosition(pos.Vect())) {
      containment_index = i;
      break;
    }
  }
  // contruct pos Point
  for (int i = 0; i < _config.active_volumes.size(); i++) {
    if (_config.active_volumes[i].ContainsPosition(pos.Vect())) {
      cryostat_index = i;
      break;
    }
  }

  // only consider contained length in the cryostat volume containing the interaction
  // setup TPC index
  std::vector<geo::BoxBoundedGeo> volumes;
  if (cryostat_index >= 0) {
    volumes = _config.tpc_volumes[cryostat_index];
    // setup the initial TPC index
    for (int i = 0; i < volumes.size(); i++) {
      if (volumes[i].ContainsPosition(pos.Vect())) {
        tpc_index = i;
        break;
      }
    }
  }
  // if we couldn't find a volume for the intial point, set not contained
  else {
    contained_in_cryo = false;
  }
  if (tpc_index < 0) {
    contained_in_tpc = false;
  }


  if (containment_index < 0) {
    is_contained = false;
  }

  // setup aa volumes too for length calc
  std::vector<geoalgo::AABox> aa_volumes;
  for (auto const &v: volumes) {
    aa_volumes.emplace_back(v.MinX(), v.MinY(), v.MinZ(), v.MaxX(), v.MaxY(), v.MaxZ());
  } 

  // Get the length and determine if any point leaves the active volume
  //
  // Use every trajectory point if possible
  if (track.NumberTrajectoryPoints() != 0) {
    // particle trajectory
    const simb::MCTrajectory &trajectory = track.Trajectory();

    for (int i = 1; i < track.NumberTrajectoryPoints(); i++) {
      TVector3 this_point = trajectory.Position(i).Vect();
      // update if track is contained
      if (contained_in_cryo) {
        contained_in_cryo = _config.active_volumes[cryostat_index].ContainsPosition(this_point);
      }
      // check if track has crossed TPC
      if (contained_in_cryo && !crosses_tpc) {
        for (int j = 0; j < volumes.size(); j++) {
          if (volumes[j].ContainsPosition(this_point) && j != tpc_index) {
            crosses_tpc = true;
            break;
          }
        }
      }
      // check if track has left tpc
      if (contained_in_tpc) {
        contained_in_tpc = volumes[tpc_index].ContainsPosition(this_point);
      }

      if (is_contained) {
        is_contained = _config.containment_volumes[containment_index].ContainsPosition(this_point);
      }
      
      // update length
      contained_length += containedLength(this_point, pos.Vect(), aa_volumes);
      pos = trajectory.Position(i);
    }
  }
  else if (_config.shakyMCTracks) {
    TVector3 end_point = track.EndPosition().Vect();

    contained_length = containedLength(track.Position().Vect(), end_point, aa_volumes);

    // track can cross TPC if start point was in TPC and Cryo
    if (contained_in_cryo && contained_in_tpc) {
      for (int i = 0; i < volumes.size(); i++) {
        if (volumes[i].ContainsPosition(end_point) && i != tpc_index) {
          crosses_tpc = true;
          break;
        }
      }
    }

    // update containment with end points
    if (contained_in_cryo) contained_in_cryo = _config.active_volumes[cryostat_index].ContainsPosition(end_point);
    if (contained_in_tpc)  contained_in_tpc = volumes[tpc_index].ContainsPosition(end_point);

  }
  else {
    std::cerr << "ERROR: unexpected bad track" << std::endl;
    assert(false);
  }

  // get the distance to the truth vertex
  // double dist_to_vertex = truth.NeutrinoSet() ? (track.Position().Vect() - truth.GetNeutrino().Nu().Position(0).Vect()).Mag(): -1;
  // TODO: set
  double dist_to_vertex = -1;

  // ret info
  numu::RecoTrack ret;

  // ret.kinetic_energy = kinetic_energy;
  ret.pdgid = pdgid;
  ret.is_muon = abs(pdgid) == 13;
  ret.length = contained_length;
  ret.costh = costh;
  ret.contained_in_cryo = contained_in_cryo;
  ret.contained_in_tpc = contained_in_tpc;
  ret.crosses_tpc = crosses_tpc;
  ret.is_contained = is_contained;
  ret.start = track.Position().Vect();
  ret.end = track.EndPosition().Vect();
  ret.momentum = track.Momentum().Vect().Mag();
  ret.energy = kinetic_energy;
  ret.dist_to_vertex = dist_to_vertex;

  return ret;
}



int NumuReco::MCTruthPrimaryTrack(const simb::MCTruth &mctruth, const std::vector<simb::MCParticle> &mcparticle_list) {
  // Get the length and determine if any point leaves the active volume
  //
  // get lepton track
  // if multiple, get the one with the highest energy
  int track_ind = -1;
  for (int i = 0; i < mcparticle_list.size(); i++) {
    if (isFromNuVertex(mctruth, mcparticle_list[i]) && abs(mcparticle_list[i].PdgCode()) == 13 && mcparticle_list[i].Process() == "primary") {
      if (track_ind == -1 || mcparticle_list[track_ind].E() < mcparticle_list[i].E()) {
        track_ind = i;
      }
    }
  }
  /*
  // if there's no lepton, look for a pi+ that can "fake" a muon
  // if there's multiple, get the one with the highest energy
  if (track_ind == -1 && mctruth.GetNeutrino().CCNC() == 1) {
    // assert(mctruth.GetNeutrino().CCNC() == 1);
    double track_contained_length = -1;
    for (int i = 0; i < mcparticle_list.size(); i++) {
      if (isFromNuVertex(mctruth, mcparticle_list[i]) && abs(mcparticle_list[i].PdgCode()) == 211 && mcparticle_list[i].Process() == "primary") {
        if (track_ind == -1 || mcparticle_list[track_ind].E() < mcparticle_list[i].E()) {
          track_ind = mcparticle_list[i].TrackId();
        }
      }
    }
  }*/
  return (track_ind == -1) ? track_ind : mcparticle_list[track_ind].TrackId();
}

int NumuReco::TrueTrackMultiplicity(const simb::MCTruth &mc_truth, const std::vector<simb::MCParticle> &mcparticle_list) {
  int multiplicity = 0;
  for (int i = 0; i < mcparticle_list.size(); i++) {
    if (isFromNuVertex(mc_truth, mcparticle_list[i]) &&  mcparticle_list[i].Process() == "primary") {
      if (PDGCharge(mcparticle_list[i].PdgCode()) > 1e-4 && mcparticle_list[i].E() - PDGMass(mcparticle_list[i].PdgCode()) > 21) {
        multiplicity ++;
      }
    }
  }
  return multiplicity;
}

std::map<size_t, numu::RecoTrack> NumuReco::MCParticleTracks(const gallery::Event &event) {
  std::map<size_t, numu::RecoTrack> ret;
  const std::vector<simb::MCParticle> &mcparticles = *event.getValidHandle<std::vector<simb::MCParticle>>(_config.MCParticleTag);
  for (unsigned i = 0; i < mcparticles.size(); i++) {
    ret[mcparticles[i].TrackId()] = std::move(MCTrackInfo(mcparticles[i]));
  }
  return ret;
}

std::vector<size_t> NumuReco::MCTruthTracks(
    std::map<size_t,  numu::RecoTrack> &true_tracks,
    const art::FindManyP<simb::MCParticle, sim::GeneratedParticleInfo> &truth_to_particles, 
    const simb::MCTruth &mc_truth, 
    int mc_truth_index) 
{

  std::vector<size_t> ret;
  std::vector<art::Ptr<simb::MCParticle>> mcparticle_list = truth_to_particles.at(mc_truth_index);
  for (unsigned i = 0; i < mcparticle_list.size(); i++) {
    true_tracks.at(mcparticle_list[i]->TrackId()).match = MCTrackMatch(mc_truth, *mcparticle_list[i]);
    ret.push_back(mcparticle_list[i]->TrackId());
  }
  return ret;
}

std::vector<numu::RecoInteraction> NumuReco::MCTruthInteractions(const gallery::Event &event, std::map<size_t, numu::RecoTrack> &true_tracks) {
  // Get truth
  gallery::Handle<std::vector<simb::MCTruth>> mctruth_handle;
  bool has_mctruth = event.getByLabel(fTruthTag, mctruth_handle);

  // also cosmics
  gallery::Handle<std::vector<simb::MCTruth>> cosmics;
  bool has_cosmics = event.getByLabel(_config.CorsikaTag, cosmics);

  auto const& mcparticle_list = \
    *event.getValidHandle<std::vector<simb::MCParticle>>(fMCParticleTag);

  std::vector<numu::RecoInteraction> ret;


  // iterate over truth interactions
  unsigned n_truth = has_mctruth ? mctruth_handle->size() :0;
  for (unsigned i = 0; i < n_truth; i++) {
    // match Genie to G4
    art::FindManyP<simb::MCParticle, sim::GeneratedParticleInfo> truth_to_particles(mctruth_handle, event, fMCParticleTag);

    auto const &truth = (*mctruth_handle)[i];
    auto neutrino = truth.GetNeutrino();
    TVector3 position = neutrino.Nu().Position().Vect();

    numu::RecoInteraction this_interaction;
    this_interaction.position = position;
    this_interaction.nu_energy = neutrino.Nu().E();
    this_interaction.slice.tracks = MCTruthTracks(true_tracks, truth_to_particles, truth, i); 
    this_interaction.multiplicity = TrueTrackMultiplicity(truth, mcparticle_list);
    // get the primary track
    this_interaction.slice.primary_track_index = MCTruthPrimaryTrack(truth, mcparticle_list);

    this_interaction.match.has_match = true;
    // get the matching info (to itself)
    if (neutrino.CCNC() == simb::kCC) {
      this_interaction.match.mode = numu::mCC;
    }
    else {
      this_interaction.match.mode = numu::mNC;
    }
    this_interaction.match.tmode = numu::tmNeutrino;
    this_interaction.match.mctruth_vertex_id = i;
    this_interaction.match.event_vertex_id = ret.size();
    this_interaction.match.mctruth_track_id = i;
    this_interaction.match.event_track_id = ret.size();
    this_interaction.match.truth_vertex_distance = 0;
    this_interaction.match.is_misreconstructed = false;

    ret.push_back(std::move(this_interaction)); 
  }
  // iterate over cosmics
  if (has_cosmics) {
    // match corsika to G4
    art::FindManyP<simb::MCParticle, sim::GeneratedParticleInfo> cosmic_to_particles(cosmics, event, fMCParticleTag);
    for (unsigned i = 0; i < cosmics->size(); i++) {
      auto const &cosmic = (*cosmics)[i];
      numu::RecoInteraction this_interaction;
      this_interaction.position = InvalidTVector3;
      this_interaction.nu_energy = -1; 
      this_interaction.slice.tracks = std::move(MCTruthTracks(true_tracks, cosmic_to_particles, cosmic, i));
      this_interaction.multiplicity = -1;
      this_interaction.slice.primary_track_index = -1;
      this_interaction.match.has_match = true;
      this_interaction.match.mode = numu::mCosmic; 
      this_interaction.match.tmode = numu::tmCosmic;
      this_interaction.match.mctruth_vertex_id = n_truth + i;
      this_interaction.match.event_vertex_id = ret.size();
      this_interaction.match.mctruth_track_id = -1;
      this_interaction.match.event_track_id = ret.size();
      this_interaction.match.truth_vertex_distance = -1;
      this_interaction.match.is_misreconstructed = false;
      ret.push_back(std::move(this_interaction));
    }
  }

  return ret;
}

double RecoTrackLength(const art::Ptr<recob::Track> &track) {
  if (track->CountValidPoints() == 0) return 0.;
  double dist = 0.;
  geo::Point_t first = track->Start();
  for (size_t i = 1; i < track->CountValidPoints(); i++) {
    geo::Point_t second = track->LocationAtPoint(i);
    dist += sqrt((second - first).Mag2());
    first = second;
  }
  return dist;
}

std::array<bool, 4> NumuReco::RecoTrackTopology(const art::Ptr<recob::Track> &track) {
  // returned info
  bool contained_in_cryo = true;
  bool contained_in_tpc = true;
  bool crosses_tpc = false; 

  bool is_contained = true;

  // start point
  geo::Point_t start = track->Start();
  // get the active volume that the start position is in
  int cryostat_index = -1;
  int tpc_index = -1;
  int containment_index = -1;
  for (int i = 0; i < _config.containment_volumes.size(); i++) {
    if (_config.containment_volumes[i].ContainsPosition(start)) {
      containment_index = i;
    }
    break;
  }
  for (int i = 0; i < _config.active_volumes.size(); i++) {
    if (_config.active_volumes[i].ContainsPosition(start)) {
      cryostat_index = i;
      break;
    }
  }
  if (containment_index < 0) {
    is_contained = false;
  }
  std::vector<geo::BoxBoundedGeo> volumes;
  if (cryostat_index >= 0) {
    volumes = _config.tpc_volumes[cryostat_index];
    for (int i = 0; i < volumes.size(); i++) {
      if (volumes[i].ContainsPosition(start)) {
        tpc_index = i;
        break;
      }
    }
  }
  else {
    contained_in_cryo = false;
  }
  if (tpc_index < 0) {
    contained_in_tpc = false;
  }

  // now check for all track points
  for (int i = 1; i < track->CountValidPoints(); i++) {
    geo::Point_t this_point = track->LocationAtPoint(i);
    if (is_contained) {
      is_contained = _config.containment_volumes[containment_index].ContainsPosition(this_point);
    }
    if (contained_in_cryo) {
      contained_in_cryo = _config.active_volumes[cryostat_index].ContainsPosition(this_point);
    }
    if (contained_in_cryo && !crosses_tpc) {
      for (int j = 0; j < volumes.size(); j++) {
        if (volumes[j].ContainsPosition(this_point) && j != tpc_index) {
          crosses_tpc = true;
          break;
        }
      }
    }
    if (contained_in_tpc) {
      contained_in_tpc = volumes[tpc_index].ContainsPosition(this_point);
    }
  }

  return {contained_in_cryo, contained_in_tpc, crosses_tpc, is_contained};
}


std::map<size_t, numu::RecoTrack> NumuReco::RecoTrackInfo() {
  std::map<size_t, numu::RecoTrack> ret;
  for (unsigned pfp_track_index = 0; pfp_track_index < _tpc_tracks.size(); pfp_track_index++) {
    const art::Ptr<recob::Track> &track = _tpc_tracks[pfp_track_index];

    // information to be saved
    numu::RecoTrack this_track;

    this_track.ID = pfp_track_index;

    // track length
    this_track.length = RecoTrackLength(track);

    // get the associated PID and Calo
    assert(_tpc_tracks_to_pid.at(pfp_track_index).size() == 3); // one per plane
    
    // sum up all the pid scores weighted by n dof
    double chi2_proton = 0.;
    double chi2_kaon = 0.;
    double chi2_muon = 0.;
    double chi2_pion = 0.;
    int n_dof = 0;
    int particle_pdg = 0;
    double min_chi2 = 0.;
    for (int i =0; i < 3; i++) {
      // invalid plane means invalid calorimetry
      if (!_tpc_tracks_to_pid.at(pfp_track_index).at(i)->PlaneID()) continue;
      const art::Ptr<anab::ParticleID> &particle_id = _tpc_tracks_to_pid.at(pfp_track_index).at(i);

      n_dof += particle_id->Ndf();
      chi2_proton += particle_id->Chi2Proton() * particle_id->Ndf();
      chi2_kaon += particle_id->Chi2Kaon() * particle_id->Ndf();
      chi2_pion += particle_id->Chi2Pion() * particle_id->Ndf();
      chi2_muon += particle_id->Chi2Muon() * particle_id->Ndf();
    }
    if (n_dof > 0) {
      /*
      chi2_proton /= n_dof;
      chi2_kaon /= n_dof;
      chi2_pion /= n_dof;
      chi2_muon /= n_dof;*/
      // min chi2 is PID
      std::vector<double> chi2s {chi2_proton, chi2_muon, chi2_kaon, chi2_pion};
      int min_ind = std::distance(chi2s.begin(), std::min_element(chi2s.begin(), chi2s.end()));
      min_chi2 = *std::min_element(chi2s.begin(), chi2s.end());
      if (min_ind == 0) {
        particle_pdg = 2212;
      }
      else if (min_ind == 1) {
        particle_pdg = 13;
      }
      else if (min_ind == 2) {
        particle_pdg = 312;
      }
      else if (min_ind == 3) {
        particle_pdg = 211;
      }
      else {
        assert(false);
      }
    }
    else {
      // No particle ID was provided -- set things to nonsense
      chi2_proton = -1;
      chi2_kaon = -1;
      chi2_muon = -1;
      chi2_pion = -1;
      min_chi2 = -1.5;
    }
    this_track.pdgid = particle_pdg;
    this_track.chi2_proton = chi2_proton;
    this_track.chi2_kaon = chi2_kaon;
    this_track.chi2_pion = chi2_pion;
    this_track.chi2_muon = chi2_muon;
    this_track.min_chi2 = min_chi2;
    this_track.pid_n_dof = n_dof;

    assert(_tpc_tracks_to_calo.at(pfp_track_index).size() == 3);
    // average and sum the deposited energies
    int n_calo = 0;
    std::vector<double> deposited_energies;
    double deposited_energy_avg = 0.;
    double deposited_energy_max = 0.;
    double deposited_energy_med = 0.;
    for (int i =0; i < 3; i++) {
      if (!_tpc_tracks_to_calo.at(pfp_track_index).at(i)->PlaneID()) continue;
      const art::Ptr<anab::Calorimetry> &calo = _tpc_tracks_to_calo.at(pfp_track_index).at(i);
      if (calo->KineticEnergy() > 1000000) {
        std::cout << "Baaaaadddd energy: " << calo->KineticEnergy() << std::endl;
        continue;
      }
      n_calo ++;
      deposited_energies.push_back(calo->KineticEnergy() / 1000.); /* MeV -> GeV */
    }
    if (n_calo > 0) {
      deposited_energy_avg = std::accumulate(deposited_energies.begin(), deposited_energies.end(), 0.) / n_calo;
      deposited_energy_max = *std::max_element(deposited_energies.begin(), deposited_energies.end());

      // 1 or 3 values -- take the middle
      if (n_calo != 2) {
        // compute the median
        std::sort(deposited_energies.begin(), deposited_energies.end());
        deposited_energy_med = deposited_energies[deposited_energies.size() / 2];
      }
      // otherwise take the average
      else if (n_calo == 2) {
        deposited_energy_med = deposited_energy_avg;
      }
      // bad
      else assert(false);
      
    }
    else {
      deposited_energy_max = -1;
      deposited_energy_avg = -1;
      deposited_energy_med = -1;
    }

    this_track.deposited_energy_max = deposited_energy_max;
    this_track.deposited_energy_avg = deposited_energy_avg;
    this_track.deposited_energy_med = deposited_energy_med;
    // calculator only has inputs for protons and muons
    int track_mom_pdg = (this_track.pdgid == 13 || this_track.pdgid == 211) ? 13 : 2212;
    this_track.range_momentum = _track_momentum_calculator->GetTrackMomentum(this_track.length, track_mom_pdg);

    this_track.range_momentum_muon = _track_momentum_calculator->GetTrackMomentum(this_track.length, 13);

    recob::MCSFitResult mcs_fit = _mcs_fitter->fitMcs(*track, this_track.pdgid);
    this_track.fwd_mcs_momentum = mcs_fit.fwdMomentum();
    this_track.fwd_mcs_momentum_err = mcs_fit.fwdMomUncertainty();
    this_track.bwd_mcs_momentum = mcs_fit.bwdMomentum();
    this_track.bwd_mcs_momentum_err = mcs_fit.bwdMomUncertainty();
    this_track.mcs_is_backward = !mcs_fit.isBestFwd();

    recob::MCSFitResult mcs_fit_muon = _mcs_fitter->fitMcs(*track, 13);
    this_track.fwd_mcs_momentum_muon = mcs_fit_muon.fwdMomentum();
    this_track.fwd_mcs_momentum_muon_err = mcs_fit_muon.fwdMomUncertainty();
    this_track.bwd_mcs_momentum_muon = mcs_fit_muon.bwdMomentum();
    this_track.bwd_mcs_momentum_muon_err = mcs_fit_muon.bwdMomUncertainty();

    // TODO: fill this
    this_track.momentum = -1;
    this_track.energy = -1;
    this_track.dist_to_vertex = -1;

    this_track.costh = track->StartDirection().Z() / sqrt( track->StartDirection().Mag2() );  

    // get track topology
    std::array<bool, 4> topology = RecoTrackTopology(track);
    this_track.contained_in_cryo = topology[0];
    this_track.contained_in_tpc = topology[1];
    this_track.crosses_tpc = topology[2];
    this_track.is_contained = topology[3];

    this_track.start = TVector3(track->Start().X(), track->Start().Y(), track->Start().Z());
    this_track.end = TVector3(track->End().X(), track->End().Y(), track->End().Z());

    // do truth matching
    this_track.match = MatchTrack2Truth(pfp_track_index);

    // if configured, apply Cosmic ID to all tracks
    if (_config.CosmicIDAllTracks) {
      ApplyCosmicID(this_track);
    }

    ret[this_track.ID] = this_track;
  }
  return ret;
}

void NumuReco::ApplyCosmicID(numu::RecoTrack &track) {
  const art::Ptr<recob::Track> &pfp_track = _tpc_tracks[track.ID];

  // do Flash Matching
  track.flash_match = FlashMatching(*pfp_track, track);
  
  // do CRT matching
  track.crt_match = CRTMatching(track, *pfp_track, _tpc_tracks_to_hits.at(track.ID));
  
  // Cosmic ID: 
  //
  // See if start or end looks like stopping
  track.stopping_chisq_start = _stopping_cosmic_alg.StoppingChiSq(pfp_track->Vertex(), _tpc_tracks_to_calo.at(track.ID));
  track.stopping_chisq_finish = _stopping_cosmic_alg.StoppingChiSq(pfp_track->End(), _tpc_tracks_to_calo.at(track.ID));
  
  // gather the pandora T0's
  // first get the particle
  unsigned particle_index = _tpc_tracks_to_particle_index.at(track.ID);
  const std::vector<art::Ptr<anab::T0>> &pandora_T0s = _tpc_particles_to_T0.at(particle_index);
  for (const art::Ptr<anab::T0> &t0: pandora_T0s) {
    track.tpc_t0s.push_back(t0->Time() * 1e-3); // convert ns -> us
  }
}


std::vector<size_t> NumuReco::RecoSliceTracks(
    const std::map<size_t, numu::RecoTrack> &tracks,
    const std::map<size_t, numu::RecoParticle> &particles) { 

  std::vector<size_t> ret;

  for (const auto &track_pair: tracks) {
    const numu::RecoTrack &track = track_pair.second;
    if (particles.count(_tpc_tracks_to_particle_index[track.ID])) {
      ret.push_back(_tpc_tracks_to_particle_index[track.ID]);
    }
  }

  return ret;
}

std::vector<numu::RecoParticle> NumuReco::RecoParticleInfo() {
  // ret
  std::vector<numu::RecoParticle> ret;

  // iterate over all of the pfparticles
  for (size_t i = 0; i < _tpc_particles.size(); i++) {
    // new reco particle
    numu::RecoParticle this_particle;
    this_particle.ID = i;

    // get the PFParticle
    const recob::PFParticle& this_pfp = *_tpc_particles.at(i);
    // get its daughters in the partcle "flow"
    this_particle.daughters.assign(_tpc_particles_to_daughters[i].begin(), _tpc_particles_to_daughters[i].end());
    // get the metadata
    const larpandoraobj::PFParticleMetadata& this_metadata = *_tpc_particles_to_metadata.at(i);
    // and the properties dict
    auto const &properties = this_metadata.GetPropertiesMap();
    // get the reco vertices
    for (const art::Ptr<recob::Vertex> vert: _tpc_particles_to_vertex.at(i)) {
      this_particle.vertices.push_back(vert->position());
    }

    // access pandora special values
    if (properties.count("IsClearCosmic")) {
      this_particle.p_is_clear_cosmic = properties.at("IsClearCosmic");
    }
    else {
      this_particle.p_is_clear_cosmic = false;
    }
    if (properties.count("NuScore")) {
      this_particle.p_nu_score = properties.at("NuScore");
    }
    else {
      this_particle.p_nu_score = -1;
    }
    if (properties.count("IsNeutrino")) {
      this_particle.p_is_neutrino = properties.at("IsNeutrino");
    }
    else {
      this_particle.p_is_neutrino = false;
    }

    ret.push_back(std::move(this_particle));
  }

  return ret;
}


int NumuReco::SelectPrimaryTrack(const std::map<size_t, numu::RecoTrack> &tracks, const numu::RecoSlice &slice) {
  // if no primary particle, then no primary track
  if (slice.primary_index < 0) return -1;

  int primary_track_index = -1;

  // try to identify the "primary" track
  // Method 0 -- "reconstruct": take the longest one with PID == 13
  // TODO: should the primary track be a daughter of the Vertex?
  // For now -- yes
  const numu::RecoParticle &neutrino = slice.particles.at(slice.primary_index);
  double max_len = -1.;
  for (size_t pfp_index: neutrino.daughters) {
    // is the daughter a track?
    if (_tpc_particles_to_track_index.count(pfp_index)) {
      unsigned track_index = _tpc_particles_to_track_index[pfp_index];
      // get it
      const numu::RecoTrack &track = tracks.at(track_index);
      if (track.is_muon) {
        if (track.length > max_len) {
          primary_track_index = track.ID;
          max_len = track.length;
        }
      }
    }
  }

  return primary_track_index;
}

std::vector<numu::RecoSlice> NumuReco::SelectSlices(const std::vector<numu::RecoSlice>& reco_slices) {
  std::vector<numu::RecoSlice> ret;

  for (auto const &slice: reco_slices) {
    if (slice.primary_index >= 0) {
      if (slice.particles.at(slice.primary_index).p_is_neutrino) {
        if (slice.primary_track_index >= 0) {
          ret.push_back(slice);
        }
      }
    }
  }

  return ret;
}

numu::TrackTruthMatch NumuReco::MatchTrack2Truth(size_t pfp_track_id) {
  // get the service
  const cheat::ParticleInventory *inventory_service = fProviderManager->GetParticleInventoryProvider();

  // get its hits
  std::vector<art::Ptr<recob::Hit>> hits = _tpc_tracks_to_hits.at(pfp_track_id);
  // this id is the same as the mcparticle ID as long as we got it from geant4
  // int id = SBNRecoUtils::TrueParticleIDFromTotalRecoHits(*fProviderManager, hits, false);
  int mcp_track_id = SBNRecoUtils::TrueParticleIDFromTotalTrueEnergy(*fProviderManager, hits, true);

  int mcparticle_index = -1;
  for (int i = 0; i < _true_particles.size(); i++) {
    if (_true_particles[i]->TrackId() == mcp_track_id) {
      mcparticle_index = i;
      break;
    }
  }

  // number returned to mean "NULL"
  // no match
  if (mcparticle_index == -1) return numu::TrackTruthMatch();

  // We got a match! Try to identify it to an origin
  art::Ptr<simb::MCTruth> truth = inventory_service->TrackIdToMCTruth_P(mcp_track_id);
  // and calculate the completion
  double completion = TrackCompletion(mcp_track_id, hits);

  numu::TrackTruthMatch ret;
  // ret.mctruth = truth.get();
  ret.has_match = true;
  ret.mctruth_has_neutrino = truth->NeutrinoSet();
  ret.mctruth_vertex = ret.mctruth_has_neutrino ? truth->GetNeutrino().Nu().Position().Vect() : TVector3(-999, -999, -999);
  ret.mctruth_origin = truth->Origin();
  ret.mctruth_ccnc = ret.mctruth_has_neutrino ? truth->GetNeutrino().CCNC() : -1;
  ret.mcparticle_id = mcp_track_id;
  ret.completion = completion;
  ret.match_pdg = _true_particles[mcparticle_index]->PdgCode(); 
  return ret;
}

std::vector<numu::RecoSlice> NumuReco::RecoSliceInfo(
    std::map<size_t, numu::RecoTrack> &reco_tracks,
    const std::vector<numu::RecoParticle> &particles) {

  std::vector<numu::RecoSlice> ret;
  // store all the particles in the slice
  for (size_t i = 0; i < _tpc_slices.size(); i++) {
    const recob::Slice &slice = *_tpc_slices[i];
    numu::RecoSlice slice_ret;
    slice_ret.primary_index = -1;
    bool primary_particle_set = false;
    // get the particles
    const std::vector<art::Ptr<recob::PFParticle>> &pfp_particles = _tpc_slices_to_particles.at(i);
    // find the primary one and store all of them
    // Find the priamry particle
    for (unsigned j = 0; j < pfp_particles.size(); j++) {
      const art::Ptr<recob::PFParticle> pfp_part = pfp_particles[j];
      if (pfp_part->IsPrimary()) {
        assert(!primary_particle_set);
        primary_particle_set = true;
        slice_ret.primary_index = _tpc_slices_to_particle_index[i][j];
      }
    }
    // match the RecoParticle's to the once that match this slice
    for (const numu::RecoParticle &part: particles) {
      if (std::find(_tpc_slices_to_particle_index[i].begin(), _tpc_slices_to_particle_index[i].end(), (unsigned)part.ID) != _tpc_slices_to_particle_index[i].end()) {
        slice_ret.particles[part.ID] = part;
      }
    } 
    
    // now get information from particles which are tracks
    slice_ret.tracks = RecoSliceTracks(reco_tracks, slice_ret.particles);
    // select the primary track
    slice_ret.primary_track_index = SelectPrimaryTrack(reco_tracks, slice_ret);
    // if we didn't do the cosmic ID for all tracks, do it for the primary
    if (slice_ret.primary_track_index >= 0 && !_config.CosmicIDAllTracks) {
      ApplyCosmicID(reco_tracks.at(slice_ret.primary_track_index));
    }

    ret.push_back(std::move(slice_ret));
  }

  return ret;
}

void NumuReco::CollectTruthInformation(const gallery::Event &ev) {
  _true_particles.clear();
  _true_particles_to_truth.clear();
  _true_particles_to_generator_info.clear();

  // mc particles
  auto const &mcparticles = ev.getValidHandle<std::vector<simb::MCParticle>>(fMCParticleTag);
  art::fill_ptr_vector(_true_particles, mcparticles);

  // MC particle (G4) to MCTruth (generator -- corsika or genie)
  art::FindManyP<simb::MCTruth, sim::GeneratedParticleInfo> particles_to_truth(mcparticles, ev, "largeant");
  for (unsigned i = 0; i < mcparticles->size(); i++) {
    _true_particles_to_truth.push_back(particles_to_truth.at(i).at(0));
    _true_particles_to_generator_info.push_back(particles_to_truth.data(i).at(0));
  }
}

void NumuReco::CollectTPCInformation(const gallery::Event &ev) {
  // clear out old stuff
  _tpc_slices.clear();
  _tpc_tracks.clear();
  _tpc_particles.clear();
  _tpc_slices_to_particles.clear();
  _tpc_slices_to_particle_index.clear();
  _tpc_tracks_to_particles.clear();
  _tpc_tracks_to_particle_index.clear();
  _tpc_particles_to_track_index.clear();
  _tpc_tracks_to_calo.clear();
  _tpc_tracks_to_pid.clear();
  _tpc_tracks_to_hits.clear();
  _tpc_particles_to_T0.clear();
  _tpc_particles_to_vertex.clear();
  _tpc_particles_to_daughters.clear();
  _tpc_particles_to_metadata.clear();

  for (const std::string &suffix: _config.TPCRecoTagSuffixes) {
    // Offset to convert a pandora particle ID to a NumuReco particle ID
    //
    // The different ID's are required because in ICARUS TPC reconstruction
    // is done per-cryostat. Thus, the pandora particle ID's may not be unique
    unsigned particle_id_offset = _tpc_particles.size();

    // get the slices
    const auto &slice_handle = ev.getValidHandle<std::vector<recob::Slice>>(_config.RecoSliceTag + suffix);
    art::fill_ptr_vector(_tpc_slices, slice_handle);

    // get the particles
    const auto &particle_handle = ev.getValidHandle<std::vector<recob::PFParticle>>(_config.PFParticleTag + suffix);
    art::fill_ptr_vector(_tpc_particles, particle_handle);
   
    // get the tracks
    const auto &track_handle = ev.getValidHandle<std::vector<recob::Track>>(_config.RecoTrackTag + suffix);
    art::fill_ptr_vector(_tpc_tracks, track_handle);

    // slice to particle
    art::FindManyP<recob::PFParticle> slice_to_particle(slice_handle, ev, _config.PFParticleTag + suffix); 
    for (unsigned i = 0; i < slice_handle->size(); i++) {
      _tpc_slices_to_particles.push_back(slice_to_particle.at(i));
      _tpc_slices_to_particle_index.emplace_back();
      for (unsigned j = 0; j < slice_to_particle.at(i).size(); j++) {
        _tpc_slices_to_particle_index[i].push_back(slice_to_particle.at(i)[j]->Self() + particle_id_offset);
        assert(_tpc_particles[_tpc_slices_to_particle_index[i][j]] == _tpc_slices_to_particles[i][j]);
      } 
    }
  
    // track to particle and particle to track
    art::FindManyP<recob::PFParticle> track_to_particle(track_handle, ev, _config.RecoTrackTag + suffix);
    for (unsigned i = 0; i < track_handle->size(); i++) {
      _tpc_tracks_to_particles.push_back(track_to_particle.at(i).at(0));
      _tpc_tracks_to_particle_index.push_back(particle_id_offset + _tpc_tracks_to_particles[i]->Self());
      _tpc_particles_to_track_index[particle_id_offset + _tpc_tracks_to_particles[i]->Self()] = i;
      // if this isn't true something went wrong
      assert(_tpc_particles[_tpc_tracks_to_particle_index[i]] == _tpc_tracks_to_particles[i]);
    }

    // track to Calo
    art::FindManyP<anab::Calorimetry> track_to_calo(track_handle, ev, _config.CaloTag + suffix);
    for (unsigned i = 0; i < track_handle->size(); i++) {
      _tpc_tracks_to_calo.push_back(track_to_calo.at(i));
    }

    // track to PID
    art::FindManyP<anab::ParticleID> track_to_pid(track_handle, ev, _config.PIDTag + suffix);
    for (unsigned i = 0; i < track_handle->size(); i++) {
      _tpc_tracks_to_pid.push_back(track_to_pid.at(i));
    }

    // track to hits
    art::FindManyP<recob::Hit> track_to_hits(track_handle, ev, _config.RecoTrackTag + suffix);
    for (unsigned i = 0; i < track_handle->size(); i++) {
      _tpc_tracks_to_hits.push_back(track_to_hits.at(i));
    }

    // particle to T0 
    art::FindManyP<anab::T0> particle_to_T0(particle_handle, ev, _config.PFParticleTag + suffix);
    for (unsigned i = 0; i < particle_handle->size(); i++) {
      _tpc_particles_to_T0.push_back(particle_to_T0.at(i));
    }

    // particle to vertex
    art::FindManyP<recob::Vertex> particle_to_vertex(particle_handle, ev, _config.PFParticleTag + suffix);
    for (unsigned i = 0; i < particle_handle->size(); i++) {
      assert(particle_to_vertex.at(i).size() == 0 || particle_to_vertex.at(i).size() == 1);
      _tpc_particles_to_vertex.push_back(particle_to_vertex.at(i));
    }

    // particle to daughters
    for (unsigned i = 0; i < particle_handle->size(); i++) {
      std::vector<unsigned> daughters;
      for (unsigned d: _tpc_particles[i]->Daughters()) {
        daughters.push_back(particle_id_offset + d); 
      }
      _tpc_particles_to_daughters.push_back(std::move(daughters));
    }

    // particle to metadata
    art::FindManyP<larpandoraobj::PFParticleMetadata> pfp_metadatas(particle_handle, ev, _config.PFParticleTag + suffix);
    for (unsigned i = 0; i < particle_handle->size(); i++) {
      _tpc_particles_to_metadata.push_back(pfp_metadatas.at(i).at(0));
    }
  }
}


void NumuReco::CollectCRTInformation(const gallery::Event &ev) {
  _crt_tracks_local.clear();
  _crt_hits_local.clear();

  // collect the CRT Tracks
  gallery::Handle<std::vector<sbnd::crt::CRTTrack>> crt_tracks_sbnd;
  bool has_crt_tracks_sbnd = ev.getByLabel(_config.CRTTrackTag, crt_tracks_sbnd);
  gallery::Handle<std::vector<icarus::crt::CRTTrack>> crt_tracks_icarus;
  // bool has_crt_tracks_icarus = ev.getByLabel(_config.CRTTrackTag, crt_tracks_icarus);
  bool has_crt_tracks_icarus = false;
  _has_crt_tracks = has_crt_tracks_sbnd || has_crt_tracks_icarus;

  if (has_crt_tracks_icarus) {
    for (unsigned i = 0; i < crt_tracks_icarus->size(); i++) {
      _crt_tracks_local.emplace_back();
      memcpy(&_crt_tracks_local[i], &crt_tracks_icarus->at(i), sizeof(sbnd::crt::CRTTrack));
    }  
  }
  _crt_tracks = (_has_crt_tracks) ?
    ( (has_crt_tracks_sbnd) ? crt_tracks_sbnd.product() : &_crt_tracks_local) : 
    NULL;

  // collect the CRT Hits
  gallery::Handle<std::vector<sbnd::crt::CRTHit>> crt_hits_sbnd;
  bool has_crt_hits_sbnd = ev.getByLabel(_config.CRTHitTag, crt_hits_sbnd);
  gallery::Handle<std::vector<icarus::crt::CRTHit>> crt_hits_icarus;
  bool has_crt_hits_icarus = ev.getByLabel(_config.CRTHitTag, crt_hits_icarus);
  _has_crt_hits = has_crt_hits_sbnd || has_crt_hits_icarus;

  // Convert ICARUS Hits to SBND Hits
  if (has_crt_hits_icarus) {
    for (unsigned i = 0; i < crt_hits_icarus->size(); i++) {
      _crt_hits_local.push_back(ICARUS2SBNDCrtHit(crt_hits_icarus->at(i)));
    }  
  }

  _crt_hits = (_has_crt_hits) ?
    ( (has_crt_hits_sbnd) ? crt_hits_sbnd.product() : &_crt_hits_local) :
    NULL;
}

std::vector<numu::CRTMatch> NumuReco::CRTMatching(
   const numu::RecoTrack &track,
   const recob::Track &pandora_track, 
   const std::vector<art::Ptr<recob::Hit>> &hits) {
//RecoTrack &track) {
  
  numu::CRTMatch match;

  // try to find a match to CRT Track -- if we have one
  int crt_id = _has_crt_tracks ? _crt_track_matchalg->GetMatchedCRTTrackId(pandora_track, hits, *_crt_tracks) : -1;

  // if this failed, see if we can get a match to a hit
  if (crt_id < 0) {
    match.has_track_match = false;
    std::pair<sbnd::crt::CRTHit, double> hit_pair = std::pair<sbnd::crt::CRTHit, double>({sbnd::crt::CRTHit(), -1});
    if (_has_crt_hits) {
      if (_config.CRTHitinOpHitRange && track.flash_match.size()) {
        const numu::FlashMatch &match = track.flash_match[0];
        double time_width = (match.match_time_width + _config.CRT2OPTimeWidth) / 2;
        std::pair<double, double> time_range; 
        time_range.first = match.match_time_first - time_width;
        time_range.second = match.match_time_first + time_width;
        int drift_dir = sbnd::TPCGeoUtil::DriftDirectionFromHits(fProviderManager->GetGeometryProvider(), hits);
        hit_pair = _crt_hit_matchalg->ClosestCRTHit(pandora_track, time_range, *_crt_hits, drift_dir); 
      }
      else {
        hit_pair = _crt_hit_matchalg->ClosestCRTHit(pandora_track, hits, *_crt_hits);
      }
    }

    double distance = hit_pair.second;
    // matching failed
    if (distance < 0) {
      return {};
    }
    // distance is too far
    if (_config.CRTHitDist > 0 && distance > _config.CRTHitDist) {
      return {};
    }
    match.has_hit_match = true;
    match.hit_distance = distance;
    match.hit = hit_pair.first;
    match.match_time = ((_config.TSMode == 0) ? match.hit.ts0_ns : match.hit.ts1_ns) / 1000.;
  }
  // Track matching worked!
  else {
    match.has_track_match = true;
    match.has_hit_match = false;
    match.track = _crt_tracks->at(crt_id);
    match.match_time = match.track.ts1_ns / 1000.;
  }

  return {match};
}

void NumuReco::CollectPMTInformation(const gallery::Event &ev) {
  // std::vector<art::Ptr<recob::OpHit>> op_hit_ptrs;
  // std::vector<recob::OpHit> op_hits;

  _op_hit_ptrs.clear();
  _op_hits_local.clear();

  return;

   // get the OpDet hits
  art::ProductID local_opdets("local_opdets");
  if (_config.MakeOpHits) {
    const std::vector<raw::OpDetWaveform> &op_waveforms = *ev.getValidHandle<std::vector<raw::OpDetWaveform>>("opdaq");
    _op_hits_local = _op_hit_maker->MakeHits(op_waveforms);
    for (unsigned i = 0; i < _op_hits_local.size(); i++) {
      _op_hit_ptrs.emplace_back(local_opdets, &_op_hits_local[i], i);
    }
  }
  else {
    gallery::Handle<std::vector<recob::OpHit>> op_hits;
    ev.getByLabel(_config.OpFlashTag, op_hits);
    art::fill_ptr_vector(_op_hit_ptrs, op_hits);
  }

}


std::vector<numu::FlashMatch> NumuReco::FlashMatching(
  const recob::Track &pandora_track,
  const numu::RecoTrack &track) {

  // don't run if back tracker is not setup
  if (fProviderManager->GetPhotonBackTrackerProvider() == NULL) {
    return {};
  }

  // see if we can do the match
  if (track.match.has_match) {
    // DumpTrueStart(ev, track.match.mcparticle_id);
    // get the correct index
    int photon_mother_id = GetPhotonMotherID(track.match.mcparticle_id);
    const std::vector<art::Ptr<recob::OpHit>> hit_matches = fProviderManager->GetPhotonBackTrackerProvider()->TrackIdToOpHits_Ps(photon_mother_id, _op_hit_ptrs);
    if (hit_matches.size() > 0) {
      // average the times weighted by PE
      double average = 0.;
      double match_pe = 0.;
      double min_time = 99999999.;
      for (const art::Ptr<recob::OpHit> op_hit: hit_matches) {
        average += op_hit->PE() * op_hit->PeakTime();
        match_pe += op_hit->PE();
        if (op_hit->PeakTime() < min_time) {
          min_time = op_hit->PeakTime();
        }
      }
      if (match_pe > 0) {
        double match_time = average / match_pe;
        numu::FlashMatch match;
        match.match_time = match_time;
        match.match_time_first = min_time;
        return {match};
      }
    }
  }
  return {};
}

numu::RecoEvent NumuReco::Reconstruct(const gallery::Event &ev, std::vector<numu::RecoInteraction> truth) {
  // setup products
  CollectCRTInformation(ev);
  CollectPMTInformation(ev);
  CollectTPCInformation(ev);
  CollectTruthInformation(ev);

  // colect PFParticle information
  std::vector<numu::RecoParticle> reco_particles = RecoParticleInfo();

  // collect track information
  std::map<size_t, numu::RecoTrack> reco_tracks = RecoTrackInfo();

  // collect Pandora slice information
  std::vector<numu::RecoSlice> reco_slices = RecoSliceInfo(reco_tracks, reco_particles);

  std::vector<numu::RecoSlice> selected_slices = SelectSlices(reco_slices);

  std::vector<numu::RecoInteraction> reco;
  for (unsigned reco_i = 0; reco_i < selected_slices.size(); reco_i++) {
    const numu::RecoParticle &neutrino = selected_slices[reco_i].particles.at(selected_slices[reco_i].primary_index);

    geo::Point_t g_pos = neutrino.vertices[0];
    TVector3 reco_position(g_pos.X(), g_pos.Y(), g_pos.Z());

    numu::RecoInteraction this_interaction;

    this_interaction.slice = selected_slices[reco_i];
    this_interaction.position = reco_position;

    this_interaction.primary_track = reco_tracks.at(this_interaction.slice.primary_track_index);
  
    // TODO: get the enrgy
    this_interaction.nu_energy = -1;

    // Track multiplicity
    this_interaction.multiplicity = this_interaction.slice.particles.at(this_interaction.slice.primary_index).daughters.size();

    // Try to match against truth events
    // Vertices
    this_interaction.match.event_vertex_id = -1;
    this_interaction.match.truth_vertex_distance = -1;
    for (int truth_i = 0; truth_i < truth.size(); truth_i++) {
      // find the closest vertex
      if (this_interaction.match.truth_vertex_distance < 0 || (truth[truth_i].position - reco_position).Mag() < this_interaction.match.truth_vertex_distance) {
        this_interaction.match.truth_vertex_distance = (truth[truth_i].position - reco_position).Mag();
        this_interaction.match.event_vertex_id = truth_i;
      }
    }

    // TODO: better way to do matching
    // only keep match if distance is less than 5cm
    if (this_interaction.match.truth_vertex_distance < 0. || this_interaction.match.truth_vertex_distance > 5.) {
      std::cout << "Vertex not matched well.\n";
    }
    else {
        std::cout << "Vertex matched to neutrino.\n";
    }

    // Match to the primary track of the event
    this_interaction.match.event_track_id = -1;
    if (this_interaction.slice.primary_track_index >= 0) {
      const numu::TrackTruthMatch &ptrack_match = reco_tracks.at(this_interaction.slice.primary_track_index).match;
      if (ptrack_match.has_match && ptrack_match.mctruth_has_neutrino) {
        for (int truth_i = 0; truth_i < truth.size(); truth_i++) {
          if ((truth[truth_i].position - ptrack_match.mctruth_vertex).Mag() < 1.) {
            this_interaction.match.event_track_id = truth_i;
          }
        }
      }
      std::cout << "Primary track index: " << this_interaction.slice.primary_track_index << std::endl;

      // return the interaction mode
      if (ptrack_match.has_match && ptrack_match.mctruth_origin == simb::Origin_t::kCosmicRay) {
        this_interaction.match.tmode = numu::tmCosmic;
        std::cout << "Track matched to cosmic.\n";
      }
      else if (ptrack_match.has_match && ptrack_match.mctruth_origin == simb::Origin_t::kBeamNeutrino) {
        this_interaction.match.tmode = numu::tmNeutrino;
        std::cout << "Track matched to neutrino.\n";
      }
      else if (ptrack_match.has_match) {
        this_interaction.match.tmode = numu::tmOther;
        std::cout << "Track matched to other.\n";
      }
      else {
        this_interaction.match.tmode = numu::tmOther;
        std::cout << "Track not matched.\n";
      }
    }
    else {
      this_interaction.match.tmode = numu::tmOther;
      std::cout << "No track to match.\n";
    }

    // Use the track to match the mode of the interaction
    if (this_interaction.slice.primary_track_index >= 0 && 
      reco_tracks.at(this_interaction.slice.primary_track_index).match.has_match) {
      const numu::TrackTruthMatch &ptrack_match = reco_tracks.at(this_interaction.slice.primary_track_index).match;
      if (ptrack_match.mctruth_origin == simb::Origin_t::kBeamNeutrino && ptrack_match.mctruth_ccnc == simb::kCC) {
        this_interaction.match.mode = numu::mCC;
      }
      else if (ptrack_match.mctruth_origin == simb::Origin_t::kBeamNeutrino && ptrack_match.mctruth_ccnc == simb::kNC) {
        this_interaction.match.mode = numu::mNC;
      }
      else if (ptrack_match.mctruth_origin == simb::Origin_t::kCosmicRay) {
        this_interaction.match.mode = numu::mCosmic;
      }
      else {
        this_interaction.match.mode = numu::mOther;
      }
    }
    // **shrug**
    else {
      this_interaction.match.mode = numu::mOther;
    }

    // bad vertex
    this_interaction.match.is_misreconstructed = this_interaction.match.tmode == numu::tmNeutrino && !this_interaction.match.has_match;

    // get the mctruth indexes from the event indexes
    if (this_interaction.match.event_vertex_id != -1) this_interaction.match.mctruth_vertex_id = truth[this_interaction.match.event_vertex_id].match.mctruth_vertex_id;
    else this_interaction.match.mctruth_vertex_id = -1;
    if (this_interaction.match.event_track_id != -1) this_interaction.match.mctruth_track_id = truth[this_interaction.match.event_track_id].match.mctruth_track_id;
    else this_interaction.match.mctruth_track_id = -1;

    // store
    reco.push_back(std::move(this_interaction));
  }

  numu::RecoEvent event;
  event.truth = std::move(truth);
  event.reco = std::move(reco);
  event.reco_tracks = std::move(reco_tracks);
 
  return std::move(event);
}


  }  // namespace SBNOsc
}  // namespace ana


DECLARE_SBN_PROCESSOR(ana::SBNOsc::NumuReco)

