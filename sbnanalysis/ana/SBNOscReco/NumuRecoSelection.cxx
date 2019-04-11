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
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"

#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/GeometryCore.h"

#include "core/Event.hh"
#include "core/Experiment.hh"
#include "NumuRecoSelection.h"
#include "../SBNOsc/Utilities.h"
#include "RecoUtils/RecoUtils.h"

#include "ubcore/LLBasicTool/GeoAlgo/GeoAABox.h"

#include "larsim/MCCheater/ParticleInventory.h"

// copied in RecoUtils here to not use art services

namespace ana {
  namespace SBNOsc {

NumuRecoSelection::NumuRecoSelection() :
  SelectionBase(),
  _event_counter(0),
  _nu_count(0),
  _selected(new std::vector<NumuRecoSelection::RecoInteraction>) {}

void NumuRecoSelection::Initialize(fhicl::ParameterSet* config) {
  if (config) {
    fhicl::ParameterSet pconfig = config->get<fhicl::ParameterSet>("NumuRecoSelection");

    std::vector<fhicl::ParameterSet> FVs = \
      pconfig.get<std::vector<fhicl::ParameterSet> >("fiducial_volumes");
    for (auto const& FV : FVs) {
      double xmin = FV.get<double>("xmin");
      double ymin = FV.get<double>("ymin");
      double zmin = FV.get<double>("zmin");
      double xmax = FV.get<double>("xmax");
      double ymax = FV.get<double>("ymax");
      double zmax = FV.get<double>("zmax");
      _config.fiducial_volumes.emplace_back(xmin, xmax, ymin, ymax, zmin, zmax);
    }

    // get the active and cryo volumes
    for (auto const &cryo: fProviderManager->GetGeometryProvider()->IterateCryostats()) {
      _config.cryostat_volumes.push_back(cryo.BoundingBox());
      geo::GeometryCore::TPC_iterator iTPC = fProviderManager->GetGeometryProvider()->begin_TPC(cryo.ID()),
                                      tend = fProviderManager->GetGeometryProvider()->end_TPC(cryo.ID());
      std::vector<geo::BoxBoundedGeo> this_tpc_volumes;
      while (iTPC != tend) {
        geo::TPCGeo const& TPC = *iTPC;
        this_tpc_volumes.push_back(TPC.ActiveBoundingBox());
        iTPC++;
      }
     _config.tpc_volumes.push_back(std::move(this_tpc_volumes));
    }

    // get the beam center
    _config.beamCenterX = pconfig.get<float>("beamCenterX", 130.);
    _config.beamCenterY = pconfig.get<float>("beamCenterY", 0.);

    _config.shakyMCTracks = pconfig.get<bool>("shakyMCTracks", false);
    _config.verbose = pconfig.get<bool>("verbose", false);

    _config.requireTrack = pconfig.get<bool>("requireTrack", false);

    _config.trackMatchContainmentCut = pconfig.get<double>("trackMatchContainmentCut", -1);

    _config.primaryTrackMethod = pconfig.get<int>("primaryTrackMethod", 0);

    _config.requireMatched = pconfig.get<bool>("requireMatched", false);
    _config.requireContained = pconfig.get<bool>("requireContained", false);

    // setup weight config
    _config.uniformWeights = pconfig.get<std::vector<std::string>>("uniformWeights", {});
    _config.constantWeight = pconfig.get<double>("constantWeight", 1.0);
    _config.cosmicWeight = pconfig.get<double>("cosmicWeight", 1.0);

    // get tag names
    _config.HitTag = config->get<std::string>("HitTag", "gaushit");
    _config.RecoTrackTag = config->get<std::string>("RecoTrackTag", "pandoraTrack");
    _config.RecoSliceTag = config->get<std::string>("RecoSliceTag", "pandora");
    _config.RecoVertexTag = config->get<std::string>("RecoVertexTag", "pandora");
    _config.PFParticleTag = config->get<std::string>("PFParticleTag", "pandora");
    _config.CorsikaTag = config->get<std::string>("CorsikaTag", "corsika");
    _config.CRTTrackTag = config->get<std::string>("CRTTrackTag", "crttrack");
    _config.MCParticleTag = config->get<std::string>("MCParticleTag", "largeant");

    std::vector<fhicl::ParameterSet> c_FVs = \
      pconfig.get<std::vector<fhicl::ParameterSet> >("containment_volumes");
    for (auto const& FV : c_FVs) {
      double xmin = FV.get<double>("xmin");
      double ymin = FV.get<double>("ymin");
      double zmin = FV.get<double>("zmin");
      double xmax = FV.get<double>("xmax");
      double ymax = FV.get<double>("ymax");
      double zmax = FV.get<double>("zmax");
      _config.containment_volumes.emplace_back(xmin, xmax, ymin, ymax, zmin, zmax);
    }
  }

  // Setup histo's for root output
  fOutputFile->cd();
  for (unsigned i = 0; i < NumuRecoSelection::nHistos; i++) {
    for (const auto mode: allModes) {
      _root_histos[i][mode].track_length = new TH1D(("track_length_" + mode2Str(mode) + "_" + histoNames[i]).c_str(), "track_length", 101, -10, 1000);
      _root_histos[i][mode].track_p = new TH1D(("track_p_" + mode2Str(mode) + "_" + histoNames[i]).c_str(), "track_p", 50, 0., 5.);
      _root_histos[i][mode].nuE = new TH1D(("nuE_" + mode2Str(mode) +"_" + histoNames[i]).c_str(), "nuE", 50, 0., 5.);
      _root_histos[i][mode].beam_center_distance = new TH1D(("beam_dist_" + mode2Str(mode) + "_" + histoNames[i]).c_str(), "beam_dist", 60, 0., 300.);
      _root_histos[i][mode].Q2 = new TH1D(("Q2_" + mode2Str(mode) + "_" + histoNames[i]).c_str(), "Q2", 50, 0., 10.);
      _root_histos[i][mode].true_contained_length = new TH1D(("tt_contained_length_" + mode2Str(mode) + "_" + histoNames[i]).c_str(), "tt_contained_length", 101, -10., 1000.);
      _root_histos[i][mode].true_track_multiplicity = new TH1D(("true_track_multiplicity_" + mode2Str(mode) + "_" + histoNames[i]).c_str(), "true_track_multiplicity", 10, 0., 10.);
      _root_histos[i][mode].crosses_tpc = new TH1D(("crosses_tpc_" + mode2Str(mode) + "_" + histoNames[i]).c_str(), "crosses_tpc", 2, -0.5, 1.5);
    }
  }

  // add branches
  fTree->Branch("reco_event", &_recoEvent);
  fTree->Branch("reco_vertices", &_selected);

  hello();
}


void NumuRecoSelection::Finalize() {
  // write out histos
  fOutputFile->cd();
  for (unsigned i = 0; i < NumuRecoSelection::nHistos; i++) {
    for (unsigned j = 0; j < NumuRecoSelection::nModes; j++) {
      _root_histos[i][j].track_length->Write();
      _root_histos[i][j].track_p->Write();
      _root_histos[i][j].nuE->Write();
      _root_histos[i][j].beam_center_distance->Write();
      _root_histos[i][j].Q2->Write();
      _root_histos[i][j].true_contained_length->Write();
      _root_histos[i][j].true_track_multiplicity->Write();
      _root_histos[i][j].crosses_tpc->Write();
    }
  }
}

std::vector<std::pair<size_t, NumuRecoSelection::TrackTruthMatch>> NumuRecoSelection::Truth2RecoTracks(const gallery::Event &ev, int mcparticle_id, const std::vector<size_t> &pfp_indices) {

  auto const& mcparticle_list = \
    *ev.getValidHandle<std::vector<simb::MCParticle>>(fMCParticleTag);

  // get the MC Particle
  const simb::MCParticle *truth_track = NULL;
  for (const simb::MCParticle &part: mcparticle_list) {
    if (mcparticle_id == part.TrackId()) {
      truth_track = &part;
      break;
    }
  }
  assert(truth_track != NULL);

  std::vector<std::pair<size_t, NumuRecoSelection::TrackTruthMatch>> ret;

  for (size_t ind: pfp_indices) {
    NumuRecoSelection::TrackTruthMatch this_match = MatchTrack2Truth(ev, ind); 
    if (this_match.mcparticle_id == mcparticle_id) {
      ret.push_back({ind, std::move(this_match)});
    }
  }

  return ret;
}

double NumuRecoSelection::TrackCompletion(int mcparticle_id, const std::vector<art::Ptr<recob::Hit>> &reco_track_hits) {
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

Event::RecoInteraction NumuRecoSelection::CoreRecoInteraction(const std::vector<Event::Interaction> &truth, const NumuRecoSelection::RecoInteraction &vertex, double weight) {
  Event::RecoInteraction ret;
  if (vertex.match.mctruth_vertex_id >= 0) {
    ret.truth_index = vertex.match.mctruth_vertex_id;
    ret.truth = truth[vertex.match.mctruth_vertex_id];
  }
  ret.reco_energy = vertex.nu_energy;
  ret.weight = weight;
  return ret;
}

bool NumuRecoSelection::ProcessEvent(const gallery::Event& ev, const std::vector<Event::Interaction> &core_truth, std::vector<Event::RecoInteraction>& reco) {
  if (_event_counter % 10 == 0) {
    std::cout << "NumuRecoSelection: Processing event " << _event_counter << " "
              << "(" << _nu_count << " neutrinos selected)"
              << std::endl;
  }
  // clear out old containers
  _selected->clear();

  bool selected = false;

  // clean up containers
  _event_counter++;

  // get the truth vertex info
  std::vector<NumuRecoSelection::RecoInteraction> truth = MCTruthInteractions(ev);

  // get the event info
  _recoEvent = Reconstruct(ev, truth);

  // select amongs the reco vertices
  for (unsigned i = 0; i < _recoEvent.reco.size(); i++) {
    const NumuRecoSelection::RecoInteraction &vertex = _recoEvent.reco[i];
    // run selection
    std::array<bool, NumuRecoSelection::nCuts> cuts = ProcessRecoCuts(_recoEvent, i);
    bool pass_selection = SelectReco(cuts);

    // compute the weight for this interaction
    double weight = 1.;
    weight *= _config.constantWeight;
    // TODO: what about cosmics?
    if (vertex.match.mctruth_vertex_id >= 0) {
      for (auto const &key: _config.uniformWeights) {
         weight *= core_truth[vertex.match.mctruth_vertex_id].weightmap.at(key)[0];
      }
    }
    if (vertex.match.mode == mCosmic) {
      weight *= _config.cosmicWeight;
    }

    if (pass_selection) {
      // store reco interaction
      _selected->push_back(vertex);
      // store sbncode reco interaction
      reco.push_back(CoreRecoInteraction(core_truth, vertex, weight));
      selected = true;
      _nu_count++;
    }

    // fill histos
    for (size_t cut_i=0; cut_i < NumuRecoSelection::nCuts; cut_i++) {
      if (cuts[cut_i]) {
        double track_length = vertex.primary_track_index >= 0 ? vertex.slice.tracks.at(vertex.primary_track_index).length: -1;
	_root_histos[1+cut_i][vertex.match.mode].track_length->Fill(track_length);
	_root_histos[1+cut_i][mAll].track_length->Fill(track_length);
        if (vertex.match.event_vertex_id >= 0 || vertex.match.event_track_id >= 0) {
          int event_id;
          int mctruth_id;
          if (vertex.match.event_vertex_id >= 0 && vertex.match.event_track_id >= 0) assert(vertex.match.event_vertex_id == vertex.match.event_track_id);
          if (vertex.match.event_vertex_id == -1)  {
            event_id = vertex.match.event_track_id;
            mctruth_id = vertex.match.mctruth_track_id;
          }
          else {
            event_id = vertex.match.event_vertex_id;
            mctruth_id = vertex.match.mctruth_vertex_id;
          }

          _root_histos[1+cut_i][vertex.match.mode].nuE->Fill(core_truth[mctruth_id].neutrino.energy);
          _root_histos[1+cut_i][mAll].nuE->Fill(core_truth[mctruth_id].neutrino.energy);

          double true_track_momentum = truth[event_id].primary_track_index >= 0 ? truth[event_id].slice.tracks.at(truth[event_id].primary_track_index).momentum.Mag() : -1;

          _root_histos[1+cut_i][vertex.match.mode].track_p->Fill(true_track_momentum);
          _root_histos[1+cut_i][mAll].track_p->Fill(true_track_momentum);

          _root_histos[1+cut_i][vertex.match.mode].Q2->Fill(core_truth[mctruth_id].neutrino.Q2);
          _root_histos[1+cut_i][mAll].Q2->Fill(core_truth[mctruth_id].neutrino.Q2);

          int crosses_tpc = truth[event_id].primary_track_index >= 0 ? truth[event_id].slice.tracks.at(truth[event_id].primary_track_index).crosses_tpc: -1;

          _root_histos[1+cut_i][vertex.match.mode].crosses_tpc->Fill(crosses_tpc);
          _root_histos[1+cut_i][mAll].crosses_tpc->Fill(crosses_tpc);

          // get the distance from the beam center
          float beam_center_distance = sqrt( (core_truth[mctruth_id].neutrino.position.X() - _config.beamCenterX) * 
                                             (core_truth[mctruth_id].neutrino.position.X() - _config.beamCenterX) +
                                             (core_truth[mctruth_id].neutrino.position.Y() - _config.beamCenterY) *
                                             (core_truth[mctruth_id].neutrino.position.Y() - _config.beamCenterY));

          _root_histos[1+cut_i][vertex.match.mode].beam_center_distance->Fill(beam_center_distance);
          _root_histos[1+cut_i][mAll].beam_center_distance->Fill(beam_center_distance);

          double length = truth[event_id].primary_track_index >= 0 ? truth[event_id].slice.tracks.at(truth[event_id].primary_track_index).length: -1;

          _root_histos[1+cut_i][vertex.match.mode].true_contained_length->Fill(length);
          _root_histos[1+cut_i][mAll].true_contained_length->Fill(length);

          _root_histos[1+cut_i][vertex.match.mode].true_track_multiplicity->Fill(truth[event_id].multiplicity);
          _root_histos[1+cut_i][mAll].true_contained_length->Fill(truth[event_id].multiplicity);

        }
      }
    }
  }
  for (int truth_i = 0; truth_i < truth.size(); truth_i++) {
    // fill in truth histos
    double track_length = truth[truth_i].primary_track_index >= 0 ? truth[truth_i].slice.tracks.at(truth[truth_i].primary_track_index).length: -1;
    _root_histos[0][truth[truth_i].match.mode].track_length->Fill(track_length);
    _root_histos[0][mAll].track_length->Fill(track_length);

    _root_histos[0][truth[truth_i].match.mode].Q2->Fill(core_truth[truth_i].neutrino.Q2);
    _root_histos[0][mAll].Q2->Fill(core_truth[truth_i].neutrino.Q2);

    int crosses_tpc = truth[truth_i].primary_track_index >= 0 ? truth[truth_i].slice.tracks.at(truth[truth_i].primary_track_index).crosses_tpc: -1;
    _root_histos[0][truth[truth_i].match.mode].crosses_tpc->Fill(crosses_tpc);
     _root_histos[0][mAll].crosses_tpc->Fill(crosses_tpc);

    _root_histos[0][truth[truth_i].match.mode].nuE->Fill(core_truth[truth_i].neutrino.energy);
    _root_histos[0][mAll].nuE->Fill(core_truth[truth_i].neutrino.energy);

    double momentum = truth[truth_i].primary_track_index >= 0 ? truth[truth_i].slice.tracks.at(truth[truth_i].primary_track_index).momentum.Mag() : -1; 
    _root_histos[0][truth[truth_i].match.mode].track_p->Fill(momentum);
    _root_histos[0][truth[truth_i].match.mode].track_p->Fill(momentum);
    // get the distance from the beam center
    float beam_center_distance = sqrt( (core_truth[truth_i].neutrino.position.X() - _config.beamCenterX) * 
                                       (core_truth[truth_i].neutrino.position.X() - _config.beamCenterX) +
                                       (core_truth[truth_i].neutrino.position.Y() - _config.beamCenterY) *
                                       (core_truth[truth_i].neutrino.position.Y() - _config.beamCenterY));
    _root_histos[0][truth[truth_i].match.mode].beam_center_distance->Fill(beam_center_distance);
    _root_histos[0][mAll].beam_center_distance->Fill(beam_center_distance);

    _root_histos[0][truth[truth_i].match.mode].true_contained_length->Fill(track_length);
    _root_histos[0][mAll].true_contained_length->Fill(track_length);

    _root_histos[0][truth[truth_i].match.mode].true_track_multiplicity->Fill(truth[truth_i].multiplicity);
    _root_histos[0][mAll].true_track_multiplicity->Fill(truth[truth_i].multiplicity);
  }

  return selected;
}

std::array<bool, NumuRecoSelection::nCuts> NumuRecoSelection::ProcessRecoCuts(const RecoEvent &event, unsigned reco_vertex_index) {
  bool is_reco = true;
  bool has_primary_track = event.reco[reco_vertex_index].primary_track_index >= 0;
  bool v_matched = event.reco[reco_vertex_index].match.event_vertex_id >= 0;

  bool t_matched = false;
  if (has_primary_track) {
    int t_mcparticle_id = event.reco[reco_vertex_index].slice.tracks.at(event.reco[reco_vertex_index].primary_track_index).match.mcparticle_id;
    for (const RecoInteraction &truth: event.truth) {
      if (truth.primary_track_index >= 0 && 
        truth.slice.tracks.at(truth.primary_track_index).match.mcparticle_id == t_mcparticle_id) {
        t_matched = true;
        break;
      } 
    }
  }
  // require completion
  t_matched = t_matched && (_config.trackMatchContainmentCut < 0 || 
    event.reco[reco_vertex_index].slice.tracks.at(event.reco[reco_vertex_index].primary_track_index).match.completion > _config.trackMatchContainmentCut);
  bool is_contained = has_primary_track && event.reco[reco_vertex_index].slice.tracks.at(event.reco[reco_vertex_index].primary_track_index).is_contained;
  return {
    is_reco,
    has_primary_track,
    v_matched,
    t_matched,
    v_matched && t_matched,
    is_contained
  };
}

bool NumuRecoSelection::SelectReco(std::array<bool, NumuRecoSelection::nCuts> &cuts) {
  return 
    cuts[0] && 
    (cuts[1] || !_config.requireTrack) &&
    (cuts[4] || !_config.requireMatched) &&
    (cuts[5] || !_config.requireContained);
}

// get information associated with track
NumuRecoSelection::RecoTrack NumuRecoSelection::MCTrackInfo(const simb::MCTruth &truth, const simb::MCParticle &track) {
  // default values
  double contained_length = 0;
  bool crosses_tpc = false;

  // If the interaction is outside the active volume, then g4 won't generate positions for the track.
  // So size == 0 => outside FV
  //
  // If size != 0, then we have to check volumes
  bool contained_in_cryo = track.NumberTrajectoryPoints() > 0;
  bool contained_in_tpc = track.NumberTrajectoryPoints() > 0;

  // other truth information
  double costh = track.Pz() / track.Momentum().Vect().Mag();
  int pdgid = track.PdgCode();
  double kinetic_energy = track.E()  - PDGMass(pdgid) / 1000.;

  // setup intial track locations
  TLorentzVector pos = track.Position();
  // get the active volume that the start position is in
  int cryostat_index = -1;
  int tpc_index = -1;
  // contruct pos Point
  for (int i = 0; i < _config.cryostat_volumes.size(); i++) {
    if (_config.cryostat_volumes[i].ContainsPosition(pos.Vect())) {
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
        contained_in_cryo = _config.cryostat_volumes[i].ContainsPosition(this_point);
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
    if (contained_in_cryo) contained_in_cryo = _config.cryostat_volumes[cryostat_index].ContainsPosition(end_point);
    if (contained_in_tpc)  contained_in_tpc = volumes[tpc_index].ContainsPosition(end_point);

  }
  else {
    std::cerr << "ERROR: unexpected bad track" << std::endl;
    assert(false);
  }

  // get the distance to the truth vertex
  double dist_to_vertex = (track.Position().Vect() - truth.GetNeutrino().Nu().Position(0).Vect()).Mag();

  // get the match info
  TrackTruthMatch track_match;
  track_match.has_match = true;
  track_match.mctruth_has_neutrino = truth.NeutrinoSet();
  track_match.mctruth_vertex = truth.NeutrinoSet() ? truth.GetNeutrino().Nu().Position().Vect() : TVector3(-999, -999, -999);
  track_match.mctruth_origin = truth.Origin();
  track_match.mctruth_ccnc = truth.NeutrinoSet() ? truth.GetNeutrino().CCNC() : -1;
  track_match.mcparticle_id = track.TrackId();
  track_match.completion = 1.;
  track_match.match_pdg = pdgid;

  // ret info
  RecoTrack ret;

  ret.kinetic_energy = kinetic_energy;
  ret.pdgid = pdgid;
  ret.is_muon = abs(pdgid) == 13;
  ret.length = contained_length;
  ret.costh = costh;
  ret.contained_in_cryo = contained_in_cryo;
  ret.contained_in_tpc = contained_in_tpc;
  ret.crosses_tpc = crosses_tpc;
  //TODO: fix
  // ret.is_contained = is_contained;
  ret.start = track.Position().Vect();
  ret.end = track.EndPosition().Vect();
  ret.momentum = track.Momentum().Vect();
  ret.dist_to_vertex = dist_to_vertex;
  ret.match = track_match;

  return ret;

}

int NumuRecoSelection::MCTruthPrimaryTrack(const simb::MCTruth &mctruth, const std::vector<simb::MCParticle> &mcparticle_list) {
  // Get the length and determine if any point leaves the active volume
  //
  // get lepton track
  // if multiple, get the one with the highest energy
  int track_ind = -1;
  for (int i = 0; i < mcparticle_list.size(); i++) {
    if (isFromNuVertex(mctruth, mcparticle_list[i]) && abs(mcparticle_list[i].PdgCode()) == 13 && mcparticle_list[i].Process() == "primary") {
      if (track_ind == -1 || mcparticle_list[track_ind].E() < mcparticle_list[i].E()) {
        track_ind = mcparticle_list[i].TrackId();
        track_ind = i;
      }
    }
  } 
  // if there's no lepton, look for a pi+ that can "fake" a muon
  // if there's multiple, get the one with the highest energy
  if (track_ind == -1) {
    assert(mctruth.GetNeutrino().CCNC() == 1);
    double track_contained_length = -1;
    for (int i = 0; i < mcparticle_list.size(); i++) {
      if (isFromNuVertex(mctruth, mcparticle_list[i]) && abs(mcparticle_list[i].PdgCode()) == 211 && mcparticle_list[i].Process() == "primary") {
        if (track_ind == -1 || mcparticle_list[track_ind].E() < mcparticle_list[i].E()) {
          track_ind = mcparticle_list[i].TrackId();
        }
      }
    }
  }

  return track_ind;
}

std::pair<std::map<size_t, NumuRecoSelection::RecoTrack>, int> NumuRecoSelection::MCTruthTracks(const simb::MCTruth &mc_truth, const std::vector<simb::MCParticle> &mcparticle_list) {

  int multiplicity = 0;
  std::map<size_t, NumuRecoSelection::RecoTrack> ret;
  for (int i = 0; i < mcparticle_list.size(); i++) {
    // TODO: better way to match MCParticles to MCNeutrino interaction???
    if (isFromNuVertex(mc_truth, mcparticle_list[i]) &&  mcparticle_list[i].Process() == "primary") {
      if (PDGCharge(mcparticle_list[i].PdgCode()) > 1e-4 && mcparticle_list[i].E() - PDGMass(mcparticle_list[i].PdgCode()) > 21) {
        multiplicity ++;
      }

      ret[mcparticle_list[i].TrackId()] = std::move(MCTrackInfo(mc_truth, mcparticle_list[i]));
    }
  }

  return {ret, multiplicity};
}

std::vector<NumuRecoSelection::RecoInteraction> NumuRecoSelection::MCTruthInteractions(const gallery::Event &event) {
  // Get truth
  auto const& mctruths = \
    *event.getValidHandle<std::vector<simb::MCTruth> >(fTruthTag);
  auto const& mcparticle_list = \
    *event.getValidHandle<std::vector<simb::MCParticle>>(fMCParticleTag);

  std::vector<NumuRecoSelection::RecoInteraction> ret;

  // iterate over truth interactions
  for (unsigned i = 0; i < mctruths.size(); i++) {
    auto const &truth = mctruths[i];
    if (!truth.NeutrinoSet()) continue;
    auto neutrino = truth.GetNeutrino();
    TVector3 position = neutrino.Nu().Position().Vect();

    if (!containedInFV(position)) continue;

    RecoInteraction this_interaction;
    this_interaction.position = position;
    this_interaction.nu_energy = neutrino.Nu().E();
    auto pair = MCTruthTracks(truth, mcparticle_list);
    this_interaction.slice.tracks = std::move(pair.first);
    this_interaction.multiplicity = pair.second;
    // get the primary track
    this_interaction.primary_track_index = MCTruthPrimaryTrack(truth, mcparticle_list);

    this_interaction.match.has_match = true;
    // get the matching info (to itself)
    if (neutrino.CCNC() == simb::kCC) {
      this_interaction.match.mode = mCC;
    }
    else {
      this_interaction.match.mode = mNC;
    }
    this_interaction.match.tmode = tmNeutrino;
    this_interaction.match.mctruth_vertex_id = i;
    this_interaction.match.event_vertex_id = ret.size();
    this_interaction.match.mctruth_track_id = i;
    this_interaction.match.event_track_id = ret.size();
    this_interaction.match.truth_vertex_distance = 0;
    this_interaction.match.is_misreconstructed = false;

    ret.push_back(std::move(this_interaction)); 
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

std::array<bool, 4> NumuRecoSelection::RecoTrackTopology(const art::Ptr<recob::Track> &track) {
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
  for (int i = 0; i < _config.cryostat_volumes.size(); i++) {
    if (_config.cryostat_volumes[i].ContainsPosition(start)) {
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
      contained_in_cryo = _config.cryostat_volumes[cryostat_index].ContainsPosition(this_point);
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

std::map<size_t, NumuRecoSelection::RecoTrack> NumuRecoSelection::RecoSliceTracks(const gallery::Event &event, const std::map<size_t, NumuRecoSelection::RecoParticle> &particles, int primary_index) {
  auto const &pfp_handle = \
    event.getValidHandle<std::vector<recob::PFParticle>>("pandora");

  // matches from pfparticles to tracks
  art::FindManyP<recob::Track> pfp_tracks(pfp_handle, event, "pandoraTrack"); 

  // matches from tracks to particle ID and calo
  // art::FindManyP<anab::Calorimetry> pfp_calo(pfp_tracks, event, "pandoraCalo");
  // art::FindManyP<anab::ParticleID> pfp_pid(pfp_tracks, event, "pandoraPID");

  std::map<size_t, NumuRecoSelection::RecoTrack> ret;

  for (auto const &this_part: particles) {
    size_t pfp_id = this_part.first;
    const NumuRecoSelection::RecoParticle &part = this_part.second;

    bool has_track = pfp_tracks.at(pfp_id).size() != 0;
    if (!has_track) continue;
    assert(pfp_tracks.at(pfp_id).size() == 1);
    // get this track
    const art::Ptr<recob::Track> &track = pfp_tracks.at(pfp_id).at(0);

    // information to be saved
    RecoTrack this_track;

    // get the associated PID and Calo
    /*
    assert(pfp_calo.at(pfp_id).size() == 1);
    assert(pfp_pid.at(pfp_id).size() == 1);
    this_track.kinetic_energy = pfp_calo.at(pfp_id).at(0)->KineticEnergy();
    this_track.pdgid = pfp_pid.at(pfp_id).at(0)->Pdg();
    */
    this_track.length = RecoTrackLength(track);
    this_track.costh = track->StartDirection().Z() / sqrt( track->StartDirection().Mag2() );  

    // get track topology
    std::array<bool, 4> topology = RecoTrackTopology(track);
    this_track.contained_in_cryo = topology[0];
    this_track.contained_in_tpc = topology[1];
    this_track.crosses_tpc = topology[2];
    this_track.is_contained = topology[3];

    this_track.start = TVector3(track->Start().X(), track->Start().Y(), track->Start().Z());
    this_track.end = TVector3(track->End().X(), track->End().Y(), track->End().Z());
    // TODO -- add these
    this_track.momentum = TVector3(-1, -1, -1);

    // get the origin of the slice
    if (primary_index >= 0) {
      this_track.dist_to_vertex = sqrt((particles.at(primary_index).vertices[0] - track->Start()).Mag2());
    }
    else {
      this_track.dist_to_vertex = -1;
    }

    this_track.match = MatchTrack2Truth(event, pfp_id);

    ret[pfp_id] = std::move(this_track);
  }

  return ret;
}

std::vector<NumuRecoSelection::RecoParticle> NumuRecoSelection::RecoParticleInfo(const gallery::Event &event) {
  // reco vertices
  auto const& reco_vertices = \
    *event.getValidHandle<std::vector<recob::Vertex>>(_config.RecoVertexTag);

  // get the PFParticles
  auto const &pfp_handle = \
    event.getValidHandle<std::vector<recob::PFParticle>>(_config.PFParticleTag);

  // use these to get all the vertices
  art::FindManyP<recob::Vertex> pfp_vertices(pfp_handle, event, "pandora");
  // and the metadata
  art::FindMany<larpandoraobj::PFParticleMetadata> pfp_metadatas(pfp_handle, event, "pandora");

  // get the CRT hits
  //gallery::ValidHandle<std::vector<sbnd::crt::CRTHit>> crt_hits = \
  //  event.getValidHandle<std::vector<sbnd::crt::CRTHit>>("crthit");

  // ret
  std::vector<NumuRecoSelection::RecoParticle> ret;

  // iterate over all of the pfparticles
  for (size_t i = 0; i < pfp_handle->size(); i++) {
    // new reco particle
    NumuRecoSelection::RecoParticle this_particle;

    // get the PFParticle
    const recob::PFParticle& this_pfp = pfp_handle->at(i);
    // get its daughters in the partcle "flow"
    this_particle.daughters.assign(this_pfp.Daughters().begin(), this_pfp.Daughters().end());
    this_particle.self = this_pfp.Self();
    // get the metadata
    const larpandoraobj::PFParticleMetadata* this_metadata = pfp_metadatas.at(i).at(0);
    // and the properties dict
    auto const &properties = this_metadata->GetPropertiesMap();
    // get the reco vertices
    for (const art::Ptr<recob::Vertex> vert: pfp_vertices.at(i)) {
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


int NumuRecoSelection::SelectPrimaryTrack(const NumuRecoSelection::RecoSlice &slice, std::vector<RecoInteraction> &truth, int truth_ind) {
  int primary_track_index = -1;

  // try to identify the "primary" track
  // Method 0 -- "reconstruct": take the longest one with PID == 13
  if (_config.primaryTrackMethod == 0) {
    // TODO: should the primary track be a daughter of the Vertex?
    // For now -- yes
    const NumuRecoSelection::RecoParticle &neutrino = slice.particles.at(slice.primary_index);
    double max_len = -1.;
    for (size_t pfp_index: neutrino.daughters) {
      // is the daughter a track?
      if (slice.tracks.find(pfp_index) != slice.tracks.end()) {  
        // get it
        const NumuRecoSelection::RecoTrack &track = slice.tracks.at(pfp_index);
        if (track.is_muon) {
          if (track.length > max_len) primary_track_index = pfp_index;
        }
      }
    }
  }
  // Method 1 -- try to get from truth primary track info
  else if (_config.primaryTrackMethod == 1) {
    double max_completion = 0.;
    if (truth_ind >= 0) {
      int match_mcparticle_id = truth[truth_ind].primary_track_index;
      assert(match_mcparticle_id >= 0);
      for (auto const &pair: slice.tracks) {
        size_t this_pfp_index = pair.first;
        const NumuRecoSelection::RecoTrack &this_track = pair.second;
        if (this_track.match.mcparticle_id == match_mcparticle_id) {
          if (this_track.match.completion > max_completion) primary_track_index = this_pfp_index;
        }
      } 
    }
  }

  return primary_track_index;
}

std::vector<NumuRecoSelection::RecoSlice> NumuRecoSelection::SelectSlices(const std::vector<NumuRecoSelection::RecoSlice>& reco_slices) {
  std::vector<NumuRecoSelection::RecoSlice> ret;

  for (auto const &slice: reco_slices) {
    if (slice.primary_index >= 0) {
      if (slice.particles.at(slice.primary_index).p_is_neutrino && containedInFV(slice.particles.at(slice.primary_index).vertices.at(0))) {
        ret.push_back(slice);
      }
    }
  }

  return ret;
}

NumuRecoSelection::TrackTruthMatch NumuRecoSelection::MatchTrack2Truth(const gallery::Event &ev, size_t pfp_index) {
  // get the service
  const cheat::ParticleInventory *inventory_service = fProviderManager->GetParticleInventoryProvider();

  // get list of mc particles
  auto const& mcparticle_list = \
    *ev.getValidHandle<std::vector<simb::MCParticle>>(fMCParticleTag);

  // get list of tracks
  auto const& recob_track_list = \
   ev.getValidHandle<std::vector<recob::Track>>("pandoraTrack");

  // get the mapping of tracks to Hits
  art::FindManyP<recob::Hit> tracks_to_hits(recob_track_list, ev, "pandoraTrack");

  // get the primary track from the daughter tracks
  auto const &pfp_handle = \
  ev.getValidHandle<std::vector<recob::PFParticle>>("pandora");
  art::FindMany<recob::Track> pfp_tracks(pfp_handle, ev, "pandoraTrack");

  // get our recob::Track object
  const recob::Track *primary_track = pfp_tracks.at(pfp_index).at(0);

  // get its hits
  std::vector<art::Ptr<recob::Hit>> hits = tracks_to_hits.at(primary_track->ID());
  // this id is the same as the mcparticle ID as long as we got it from geant4
  // int id = SBNRecoUtils::TrueParticleIDFromTotalRecoHits(*fProviderManager, hits, false);
  int mcp_track_id = SBNRecoUtils::TrueParticleIDFromTotalTrueEnergy(*fProviderManager, hits, false);

  int mcparticle_index = -1;
  for (int i = 0; i < mcparticle_list.size(); i++) {
    if (mcparticle_list[i].TrackId() == mcp_track_id) {
      mcparticle_index = i;
      break;
    }
  }

  // number returned to mean "NULL"
  // no match
  if (mcparticle_index == -1) return NumuRecoSelection::TrackTruthMatch();

  // We got a match! Try to identify it to an origin
  art::Ptr<simb::MCTruth> truth = inventory_service->TrackIdToMCTruth_P(mcp_track_id);
  // and calculate the completion
  double completion = NumuRecoSelection::TrackCompletion(mcp_track_id, hits);

  TrackTruthMatch ret;
  // ret.mctruth = truth.get();
  ret.has_match = true;
  ret.mctruth_has_neutrino = truth->NeutrinoSet();
  ret.mctruth_vertex = ret.mctruth_has_neutrino ? truth->GetNeutrino().Nu().Position().Vect() : TVector3(-999, -999, -999);
  ret.mctruth_origin = truth->Origin();
  ret.mctruth_ccnc = truth->GetNeutrino().CCNC();
  ret.mcparticle_id = mcp_track_id;
  ret.completion = completion;
  ret.match_pdg = mcparticle_list[mcparticle_index].PdgCode();
  return ret;
}

std::vector<NumuRecoSelection::RecoSlice> NumuRecoSelection::RecoSliceInfo(const gallery::Event &ev, const std::vector<NumuRecoSelection::RecoParticle> &particles) {
  // get the handle to the reco slices
  auto const &reco_slices_handle = \
    ev.getValidHandle<std::vector<recob::Slice>>(_config.RecoSliceTag);

  // get the association between these and particles
  art::FindManyP<recob::PFParticle> slice_to_particles(reco_slices_handle, ev, "pandora");

  std::vector<NumuRecoSelection::RecoSlice> ret;
  // store all the particles in the slice
  for (size_t i = 0; i < reco_slices_handle->size(); i++) {
    const recob::Slice &slice = (*reco_slices_handle)[i];
    NumuRecoSelection::RecoSlice slice_ret;
    slice_ret.primary_index = -1;
    bool primary_particle_set = false;
    // get the particles
    const std::vector<art::Ptr<recob::PFParticle>> &pfp_particles = slice_to_particles.at(i);
    // find the primary one and store all of them
    for (const art::Ptr<recob::PFParticle> &pfp_part: pfp_particles) {
      if (pfp_part->IsPrimary()) {
        assert(!primary_particle_set);
        primary_particle_set = true;
        slice_ret.primary_index = pfp_part->Self();
      }
      for (const NumuRecoSelection::RecoParticle &part: particles) {
        if (part.self == pfp_part->Self()) {
          slice_ret.particles[part.self] = part;
          break;
        }
      }
    }
    // now get information from particles which are tracks
    slice_ret.tracks = RecoSliceTracks(ev, slice_ret.particles, slice_ret.primary_index);

    // assert(primary_particle_set);
    ret.push_back(std::move(slice_ret));
  }

  // Print info on the slice:
  /*
  for (auto const &slice: ret) {
    std::cout << "Primary particle: " << slice.primary_index << std::endl;
    if (slice.primary_index >= 0) {
    std::cout << "Is neutrino: " << slice.particles.at(slice.primary_index).p_is_neutrino << std::endl;
    std::cout << "Vertices: " << slice.particles.at(slice.primary_index).vertices.size() << std::endl;
    std::cout << "V0 at: " << slice.particles.at(slice.primary_index).vertices[0].X() << " " <<  slice.particles.at(slice.primary_index).vertices[0].Y() << " " <<  slice.particles.at(slice.primary_index).vertices[0].Z() << std::endl;
    std::cout << "Particle Flow:" << std::endl;
    std::queue<std::pair<size_t, size_t>> flow;
    flow.push({slice.primary_index, slice.primary_index});
    while (flow.size() != 0) {
      std::pair<size_t, size_t> this_pair = flow.front(); flow.pop();
      size_t this_parent = this_pair.first;
      size_t this_vert = this_pair.second;
      std::cout << this_parent << " " << this_vert << std::endl;
      for (size_t c: slice.particles.at(this_vert).daughters) flow.push({this_vert, c});
    }
    std::cout << std::endl << std::endl;
    }
    else {
      for (auto const &c: slice.particles) std::cout << c.first << std::endl;

    }
  }

  for (auto const &part: particles) {
    std::cout << "Part: " << part.self << std::endl;
  }*/


  return ret;
}


NumuRecoSelection::RecoEvent NumuRecoSelection::Reconstruct(const gallery::Event &ev, std::vector<NumuRecoSelection::RecoInteraction> truth) {
  std::cout << std::endl << "New Event: " << _event_counter << std::endl;
  std::vector<NumuRecoSelection::RecoParticle> reco_particles = RecoParticleInfo(ev);

  std::vector<NumuRecoSelection::RecoSlice> reco_slices = RecoSliceInfo(ev, reco_particles);

  std::vector<NumuRecoSelection::RecoSlice> selected_slices = SelectSlices(reco_slices);
 
  std::vector<NumuRecoSelection::RecoInteraction> reco;
  for (unsigned reco_i = 0; reco_i < selected_slices.size(); reco_i++) {
    const RecoParticle &neutrino = selected_slices[reco_i].particles.at(selected_slices[reco_i].primary_index);

    geo::Point_t g_pos = neutrino.vertices[0];
    TVector3 reco_position(g_pos.X(), g_pos.Y(), g_pos.Z());

    NumuRecoSelection::RecoInteraction this_interaction;

    this_interaction.slice = selected_slices[reco_i];
    this_interaction.position = reco_position;

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
      // TODO: better way to do matching
      // only keep match if distance is less than 5cm
      if (this_interaction.match.truth_vertex_distance > 5.) {
        this_interaction.match.event_vertex_id = -1;
        std::cout << "Vertex not matched.\n";
      }
      else 
          std::cout << "Vertex matched to neutrino.\n";
    }

    // select the primary track
    this_interaction.primary_track_index = SelectPrimaryTrack(this_interaction.slice, truth, this_interaction.match.event_vertex_id);
    // TODO: get the enrgy
    this_interaction.nu_energy = -1;

    // Track multiplicity
    this_interaction.multiplicity = this_interaction.slice.particles.at(this_interaction.slice.primary_index).daughters.size();

    // determine the index and the mode
    this_interaction.match.event_track_id = -1;
    if (this_interaction.primary_track_index >= 0) {
      const NumuRecoSelection::TrackTruthMatch &ptrack_match = this_interaction.slice.tracks.at(this_interaction.primary_track_index).match;
      if (ptrack_match.has_match && ptrack_match.mctruth_has_neutrino) {
        for (int truth_i = 0; truth_i < truth.size(); truth_i++) {
          if ((truth[truth_i].position - ptrack_match.mctruth_vertex).Mag() < 1.) {
            this_interaction.match.event_track_id = truth_i;
          }
        }
      }

      // return the interaction mode
      if (ptrack_match.has_match && ptrack_match.mctruth_origin == simb::Origin_t::kCosmicRay) {
        this_interaction.match.tmode = tmCosmic;
        std::cout << "Track matched to cosmic.\n";
      }
      else if (ptrack_match.has_match && ptrack_match.mctruth_origin == simb::Origin_t::kBeamNeutrino) {
        this_interaction.match.tmode = tmNeutrino;
        std::cout << "Track matched to neutrino.\n";
      }
      else {
        this_interaction.match.tmode = tmOther;
        std::cout << "Track not matched.\n";
      }
    }
    else {
      this_interaction.match.tmode = tmOther;
      std::cout << "No track to match.\n";
    }

    // determine the mode of this interaction
    //
    // Best option: compare the vertex
    if (this_interaction.match.event_vertex_id != -1) {
      this_interaction.match.mode = truth[this_interaction.match.event_vertex_id].match.mode;
    }
    // next -- try the track
    else if (this_interaction.primary_track_index && 
      this_interaction.slice.tracks.at(this_interaction.primary_track_index).match.has_match) {
      const NumuRecoSelection::TrackTruthMatch &ptrack_match = this_interaction.slice.tracks.at(this_interaction.primary_track_index).match;
      if (ptrack_match.mctruth_origin == simb::Origin_t::kBeamNeutrino && ptrack_match.mctruth_ccnc == simb::kCC) {
        this_interaction.match.mode = mCC;
      }
      else if (ptrack_match.mctruth_origin == simb::Origin_t::kBeamNeutrino && ptrack_match.mctruth_ccnc == simb::kNC) {
        this_interaction.match.mode = mNC;
      }
      else if (ptrack_match.mctruth_origin == simb::Origin_t::kCosmicRay) {
        this_interaction.match.mode = mCosmic;
      }
      else {
        this_interaction.match.mode = mOther;
      }
    }
    // **shrug**
    else {
      this_interaction.match.mode = mOther;
    }

    // bad vertex
    this_interaction.match.is_misreconstructed = this_interaction.match.tmode == tmNeutrino && !this_interaction.match.has_match;

    // get the mcturth indexes from the event indexes
    if (this_interaction.match.event_vertex_id != -1) this_interaction.match.mctruth_vertex_id = truth[this_interaction.match.event_vertex_id].match.mctruth_vertex_id;
    else this_interaction.match.mctruth_vertex_id = -1;
    if (this_interaction.match.event_track_id != -1) this_interaction.match.mctruth_track_id = truth[this_interaction.match.event_track_id].match.mctruth_track_id;
    else this_interaction.match.mctruth_track_id = -1;

    // store
    reco.push_back(std::move(this_interaction));
  }

  RecoEvent event;
  event.truth = std::move(truth);
  event.reco = std::move(reco);
  return std::move(event);
}

bool NumuRecoSelection::containedInFV(const geo::Point_t &v) {
  for (auto const& FV: _config.fiducial_volumes) {
    if (FV.ContainsPosition(v)) return true;
  }
  return false;
}

bool NumuRecoSelection::containedInFV(const TVector3 &v) {
  for (auto const& FV: _config.fiducial_volumes) {
    if (FV.ContainsPosition(v)) return true;
  }
  return false;
}

  }  // namespace SBNOsc
}  // namespace ana


DECLARE_SBN_PROCESSOR(ana::SBNOsc::NumuRecoSelection)

