#include <list>
#include <string>
#include <iostream>
#include <cmath>
#include <vector>

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
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"

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
  _selected(new std::vector<NumuRecoSelection::RecoVertex>) {}

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

    _config.requireMatched = pconfig.get<bool>("requireMatched", false);
    _config.requireContained = pconfig.get<bool>("requireContained", false);

    // setup weight config
    _config.uniformWeights = pconfig.get<std::vector<std::string>>("uniformWeights", {});
    _config.constantWeight = pconfig.get<double>("constantWeight", 1.0);
    _config.cosmicWeight = pconfig.get<double>("cosmicWeight", 1.0);

    // get tag names
    _config.HitTag = config->get<std::string>("HitTag", "gaushit");
    _config.RecoTrackTag = config->get<std::string>("RecoTrackTag", "pandoraTrack");
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
      _root_histos[i][mode].Q2 = new TH1D(("Q2_" + mode2Str(mode) + "_" + histoNames[i]).c_str(), "Q2", 50, 0., 25.);
      _root_histos[i][mode].true_contained_length = new TH1D(("tt_contained_length_" + mode2Str(mode) + "_" + histoNames[i]).c_str(), "tt_contained_length", 101, -10., 1000.);
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
    }
  }
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

Event::RecoInteraction NumuRecoSelection::CoreRecoInteraction(const std::vector<Event::Interaction> &truth, const NumuRecoSelection::RecoVertex &vertex, double weight) {
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
  std::vector<NumuRecoSelection::RecoVertex> truth = MCTruthRecoVertexInfo(ev); 
  // get the event info
  _recoEvent = Reconstruct(ev, truth);

  // select amongs the reco vertices
  for (unsigned i = 0; i < _recoEvent.reco.size(); i++) {
    const NumuRecoSelection::RecoVertex &vertex = _recoEvent.reco[i];
    // run selection
    std::array<bool, NumuRecoSelection::nCuts> cuts = ProcessRecoCuts(_recoEvent, i);
    bool pass_selection = SelectReco(cuts);

    // compute the weight for this interaction
    double weight = 1.;
    weight *= _config.constantWeight;
    // TODO: what about cosmics?
    if (vertex.match.mctruth_vertex_id >= 0) {
      for (auto const &key: _config.uniformWeights) {
         weight *= core_truth[vertex.match.mctruth_vertex_id].weights.at(key)[0];
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
	_root_histos[1+cut_i][vertex.match.mode].track_length->Fill(vertex.track.length);
	_root_histos[1+cut_i][mAll].track_length->Fill(vertex.track.length);
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

          _root_histos[1+cut_i][vertex.match.mode].track_p->Fill(truth[event_id].track.momentum.Mag());
          _root_histos[1+cut_i][mAll].track_p->Fill(truth[event_id].track.momentum.Mag());

          _root_histos[1+cut_i][vertex.match.mode].Q2->Fill(core_truth[mctruth_id].neutrino.Q2);
          _root_histos[1+cut_i][mAll].Q2->Fill(core_truth[mctruth_id].neutrino.Q2);

          // get the distance from the beam center
          float beam_center_distance = sqrt( (core_truth[mctruth_id].neutrino.position.X() - _config.beamCenterX) * 
                                             (core_truth[mctruth_id].neutrino.position.X() - _config.beamCenterX) +
                                             (core_truth[mctruth_id].neutrino.position.Y() - _config.beamCenterY) *
                                             (core_truth[mctruth_id].neutrino.position.Y() - _config.beamCenterY));

          _root_histos[1+cut_i][vertex.match.mode].beam_center_distance->Fill(beam_center_distance);
          _root_histos[1+cut_i][mAll].beam_center_distance->Fill(beam_center_distance);

          _root_histos[1+cut_i][vertex.match.mode].true_contained_length->Fill(truth[event_id].track.length);
          _root_histos[1+cut_i][mAll].true_contained_length->Fill(truth[event_id].track.length);

        }
      }
    }
  }
  for (int truth_i = 0; truth_i < truth.size(); truth_i++) {
    // fill in truth histos
    _root_histos[0][truth[truth_i].match.mode].track_length->Fill(truth[truth_i].track.length);
    _root_histos[0][mAll].track_length->Fill(truth[truth_i].track.length);

    _root_histos[0][truth[truth_i].match.mode].nuE->Fill(core_truth[truth_i].neutrino.energy);
    _root_histos[0][mAll].nuE->Fill(core_truth[truth_i].neutrino.energy);

    _root_histos[0][truth[truth_i].match.mode].track_p->Fill(truth[truth_i].track.momentum.Mag());
    _root_histos[0][mAll].track_p->Fill(truth[truth_i].track.momentum.Mag());
    // get the distance from the beam center
    float beam_center_distance = sqrt( (core_truth[truth_i].neutrino.position.X() - _config.beamCenterX) * 
                                       (core_truth[truth_i].neutrino.position.X() - _config.beamCenterX) +
                                       (core_truth[truth_i].neutrino.position.Y() - _config.beamCenterY) *
                                       (core_truth[truth_i].neutrino.position.Y() - _config.beamCenterY));
    _root_histos[0][truth[truth_i].match.mode].beam_center_distance->Fill(beam_center_distance);
    _root_histos[0][mAll].beam_center_distance->Fill(beam_center_distance);

    _root_histos[0][truth[truth_i].match.mode].true_contained_length->Fill(truth[truth_i].track.length);
    _root_histos[0][mAll].true_contained_length->Fill(truth[truth_i].track.length);
  }

  return selected;
}

std::array<bool, NumuRecoSelection::nCuts> NumuRecoSelection::ProcessRecoCuts(const RecoEvent &event, unsigned reco_vertex_index) {
  bool is_reco = true;
  bool has_primary_track = event.reco[reco_vertex_index].track.length > 0;
  bool v_matched = event.reco[reco_vertex_index].match.event_vertex_id >= 0;
  bool t_matched = event.reco[reco_vertex_index].match.event_track_id >= 0 &&
        (event.reco[reco_vertex_index].match.track_completion > _config.trackMatchContainmentCut || _config.trackMatchContainmentCut < 0.);
  bool is_contained = event.reco[reco_vertex_index].track.is_contained;
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
    (cuts[2] || !_config.requireMatched) &&
    (cuts[3] || !_config.requireMatched) &&
    (cuts[4] || !_config.requireContained);
}

// get information associated with track
NumuRecoSelection::TrackInfo NumuRecoSelection::MCTrackInfo(const simb::MCTruth &truth, const sim::MCTrack &track) {
  // default values
  double contained_length = 0;
  bool crosses_tpc = false;

  // If the interaction is outside the active volume, then g4 won't generate positions for the track.
  // So size == 0 => outside FV
  //
  // If size != 0, then we have to check volumes
  bool contained_in_cryo = track.size() > 0;
  bool contained_in_tpc = track.size() > 0;

  // other truth information
  double costh = track.Start().Momentum().Vect().Z() / track.Start().Momentum().Vect().Mag();
  int pdgid = track.PdgCode();
  double energy = track.Start().E();

  // setup intial track locations
  TLorentzVector pos = track.Start().Position();
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
  if (track.size() != 0) {
    for (int i = 1; i < track.size(); i++) {
      TVector3 this_point = track[i].Position().Vect();
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
      contained_length += containedLength(track[i].Position().Vect(), pos.Vect(), aa_volumes);
      pos = track[i].Position();
    }
  }
  else if (_config.shakyMCTracks) {
    TVector3 end_point = track.End().Position().Vect();

    contained_length = containedLength(track.Start().Position().Vect(), end_point, aa_volumes);

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

  return {
    contained_length,
    energy,
    costh,
    pdgid,
    contained_in_cryo,
    contained_in_tpc,
    crosses_tpc,
    contained_in_cryo,
    track.Start().Momentum().Vect(),
    track.End().Momentum().Vect()
  };
}

int NumuRecoSelection::GetNeutrinoTrack(const gallery::Event &ev, const simb::MCTruth &mctruth) {
  // get handle to tracks and showers
  auto const& mctrack_list = \
    *ev.getValidHandle<std::vector<sim::MCTrack> >(fMCTrackTag);

  // Get the length and determine if any point leaves the active volume
  //
  // get lepton track
  // if multiple, get the one with the highest energy
  int track_ind = -1;
  for (int i = 0; i < mctrack_list.size(); i++) {
    if (isFromNuVertex(mctruth, mctrack_list[i]) && abs(mctrack_list[i].PdgCode()) == 13 && mctrack_list[i].Process() == "primary") {
      if (track_ind == -1 || mctrack_list[track_ind].Start().E() < mctrack_list[i].Start().E()) {
        track_ind = i;
      }
    }
  } 
  // if there's no lepton, look for a pi+ that can "fake" a muon
  // if there's multiple, get the one with the highest energy
  if (track_ind == -1) {
    double track_contained_length = -1;
    for (int i = 0; i < mctrack_list.size(); i++) {
      if (isFromNuVertex(mctruth, mctrack_list[i]) && abs(mctrack_list[i].PdgCode()) == 211 && mctrack_list[i].Process() == "primary") {
        if (track_ind == -1 || mctrack_list[track_ind].Start().E() < mctrack_list[i].Start().E()) {
          track_ind = i;
        }
      }
    }
  }
  return track_ind;
}

NumuRecoSelection::TrackInfo NumuRecoSelection::MCTruthTrackInfo(const gallery::Event &event, const simb::MCTruth &mc_truth) {
  // get track info
  auto const& mctrack_list = \
    *event.getValidHandle<std::vector<sim::MCTrack> >(fMCTrackTag);

  int track_ind = GetNeutrinoTrack(event, mc_truth);
  // if there's no track, return nonsense
  if (track_ind == -1) {
    return {-1, -1, -1, -1, false, false};
  }
  else {
    return MCTrackInfo(mc_truth, mctrack_list[track_ind]);
  }
}

std::vector<NumuRecoSelection::RecoVertex> NumuRecoSelection::MCTruthRecoVertexInfo(const gallery::Event &event) {
  // Get truth
  auto const& mctruths = \
    *event.getValidHandle<std::vector<simb::MCTruth> >(fTruthTag);

  std::vector<NumuRecoSelection::RecoVertex> ret;

  // iterate over truth interactions
  for (unsigned i = 0; i < mctruths.size(); i++) {
    auto const &truth = mctruths[i];
    if (!truth.NeutrinoSet()) continue;
    auto neutrino = truth.GetNeutrino();
    TVector3 position = neutrino.Nu().Position().Vect();

    if (!containedInFV(position)) continue;

    RecoVertex this_interaction;
    this_interaction.position = position;
    this_interaction.nu_energy = neutrino.Nu().E();
    this_interaction.track = MCTruthTrackInfo(event, truth);

    // get the matching info (to itself)
    if (neutrino.CCNC() == simb::kCC) {
      this_interaction.match.mode = mCC;
    }
    else {
      this_interaction.match.mode = mNC;
    }
    this_interaction.match.tmode = tmNeutrino;
    this_interaction.match.event_vertex_id = ret.size();
    this_interaction.match.mctruth_vertex_id = i;
    this_interaction.match.mctruth_track_id = i;
    this_interaction.match.truth_vertex_distance = 0;
    this_interaction.match.is_misreconstructed = false;

    ret.push_back(std::move(this_interaction)); 
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
  art::FindMany<recob::Vertex> pfp_vertices(pfp_handle, event, "pandora");
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
    // get the metadata
    const larpandoraobj::PFParticleMetadata* this_metadata = pfp_metadatas.at(i).at(0);
    // and the properties dict
    auto const &properties = this_metadata->GetPropertiesMap();
    // get the reco vertices
    this_particle.vertices = pfp_vertices.at(i);

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

double RecoTrackLength(const recob::Track *track) {
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

std::array<bool, 4> NumuRecoSelection::RecoTrackTopology(const recob::Track *track) {
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


std::pair<NumuRecoSelection::TrackInfo, int> NumuRecoSelection::RecoTrackInfo(const gallery::Event &ev, const NumuRecoSelection::RecoParticle& vertex) {
  // get the primary track from the daughter tracks
  auto const &pfp_handle = \
    ev.getValidHandle<std::vector<recob::PFParticle>>("pandora");
  art::FindMany<recob::Track> pfp_tracks(pfp_handle, ev, "pandoraTrack");

  int primary_track_index = -1;
  double primary_track_length = 0;

  // try to identify the "primary" track
  // For now -- just take the longest one
  for (size_t i = 0; i < vertex.daughters.size(); i++) {
    size_t daughter = vertex.daughters[i];
    // get the track we care about
    const std::vector<const recob::Track *> &daughter_tracks = pfp_tracks.at(daughter);

    if (daughter_tracks.size() > 0) {
      double this_length = RecoTrackLength(daughter_tracks.at(0));
      if (this_length > 1e-4 && (primary_track_index == -1 || this_length > primary_track_length)) {
        primary_track_index = i;
        primary_track_length = this_length;
      } 
    }
  }
  // if we have identified a primary track, get its info
  if (primary_track_index != -1) {
    const recob::Track *primary_track = pfp_tracks.at(vertex.daughters[primary_track_index]).at(0);
    double length = primary_track_length;
    double costh = primary_track->StartDirection().Z() / sqrt( primary_track->StartDirection().Mag2() );  
    TVector3 start(primary_track->Start().X(), primary_track->Start().Y(), primary_track->Start().Z());
    TVector3 end(primary_track->End().X(), primary_track->End().Y(), primary_track->End().Z());

    // TODO -- add these
    double energy = -1;
    int pdgid = -1;

    // get track topology
    std::array<bool, 4> topology = RecoTrackTopology(primary_track);
    bool contained_in_cryo = topology[0];
    bool contained_in_tpc = topology[1];
    bool crosses_tpc = topology[2];
    bool is_contained = topology[3];

    // construct and return track information
    return {{
      length,
      energy,
      costh,
      pdgid,
      contained_in_cryo,
      contained_in_tpc,
      crosses_tpc,
      is_contained,
      start,
      end,
    }, primary_track_index};
  }
  // otherwise, retrun nonsense
  else {
    return {{-1, -1, -1, -1, false, false, false, -1}, -1};
  }
}

std::vector<NumuRecoSelection::RecoParticle> NumuRecoSelection::SelectVertices(const std::vector<NumuRecoSelection::RecoParticle>& reco_particles) {
  std::vector<NumuRecoSelection::RecoParticle> ret;
  for (auto const &particle: reco_particles) {
     if (particle.p_is_neutrino && containedInFV(particle.vertices[0]->position())) {
       ret.push_back(particle);
     }
  }
  return ret;
}

NumuRecoSelection::TrackTruthMatch NumuRecoSelection::MatchTrack2Truth(const gallery::Event &ev, NumuRecoSelection::RecoParticle neutrino, int track_id) {
  // if no primary track, no match
  if (track_id == -1) return NumuRecoSelection::TrackTruthMatch();

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
  size_t pfp_reco_track_index = neutrino.daughters[track_id];
  const recob::Track *primary_track = pfp_tracks.at(pfp_reco_track_index).at(0);

  // get its hits
  std::vector<art::Ptr<recob::Hit>> hits = tracks_to_hits.at(primary_track->ID());
  // this id is the same as the mcparticle ID as long as we got it from geant4
  int id = SBNRecoUtils::TrueParticleIDFromTotalRecoHits(*fProviderManager, hits, false);

  // number returned to mean "NULL"
  // no match
  if (id == -99999 || id == 0) return NumuRecoSelection::TrackTruthMatch();

  // We got a match! Try to identify it to an origin
  art::Ptr<simb::MCTruth> truth = inventory_service->TrackIdToMCTruth_P(id);
  // and calculate the completion
  double completion = NumuRecoSelection::TrackCompletion(id, hits);

  TrackTruthMatch ret;
  ret.mctruth = truth.get();
  ret.mcparticle_id = id;
  ret.completion = completion;
  return ret;
}

NumuRecoSelection::RecoEvent NumuRecoSelection::Reconstruct(const gallery::Event &ev, std::vector<NumuRecoSelection::RecoVertex> truth) {
  std::cout << std::endl << "New Event" << std::endl;
  std::vector<NumuRecoSelection::RecoParticle> reco_particles = RecoParticleInfo(ev);
  std::vector<NumuRecoSelection::RecoParticle> selected_particles = SelectVertices(reco_particles);
 
  std::vector<RecoVertex> reco;
  for (unsigned reco_i = 0; reco_i < selected_particles.size(); reco_i++) {
    const RecoParticle &particle = selected_particles[reco_i];
    assert(particle.vertices.size() == 1);

    geo::Point_t g_pos = particle.vertices[0]->position();
    TVector3 reco_position(g_pos.X(), g_pos.Y(), g_pos.Z());

    NumuRecoSelection::RecoVertex this_vertex;
    this_vertex.position = reco_position;

    // get the reco track info
    auto ret_pair = RecoTrackInfo(ev, particle);
    this_vertex.track = ret_pair.first;
    int daughter_id = ret_pair.second;
    // TODO: get the enrgy
    this_vertex.nu_energy = -1;

    // Try to match against truth events
    // Vertices
    this_vertex.match.event_vertex_id = -1;
    this_vertex.match.truth_vertex_distance = -1;
    for (int truth_i = 0; truth_i < truth.size(); truth_i++) {
      // find the closest vertex
      if (this_vertex.match.truth_vertex_distance < 0 || (truth[truth_i].position - reco_position).Mag() < this_vertex.match.truth_vertex_distance) {
        this_vertex.match.truth_vertex_distance = (truth[truth_i].position - reco_position).Mag();
        this_vertex.match.event_vertex_id = truth_i;
      }
      // TODO: better way to do matching
      // only keep match if distance is less than 5cm
      if (this_vertex.match.truth_vertex_distance > 5.) {
        this_vertex.match.event_vertex_id = -1;
        std::cout << "Vertex not matched.\n";
      }
      else 
          std::cout << "Vertex matched to neutrino.\n";
    }
    // Match against tracks
    NumuRecoSelection::TrackTruthMatch track_truth_match = MatchTrack2Truth(ev, particle, daughter_id);
    const simb::MCTruth *track_truth = track_truth_match.mctruth;

    // determine the index and the mode
    this_vertex.match.event_track_id = -1;
    if (track_truth != NULL && track_truth->NeutrinoSet()) {
      for (int truth_i = 0; truth_i < truth.size(); truth_i++) {
        if ((truth[truth_i].position - track_truth->GetNeutrino().Nu().Position(0).Vect()).Mag() < 1.) {
          this_vertex.match.event_track_id = truth_i;
        }
      }
    }
    // return the interaction mode
    if (track_truth != NULL && track_truth->Origin() ==  simb::Origin_t::kCosmicRay) {
      this_vertex.match.tmode = tmCosmic;
      std::cout << "Track matched to cosmic.\n";
    }
    else if (track_truth != NULL && track_truth->Origin() == simb::Origin_t::kBeamNeutrino) {
      this_vertex.match.tmode = tmNeutrino;
      std::cout << "Track matched to neutrino.\n";
    }
    else {
      this_vertex.match.tmode = tmOther;
      std::cout << "Track not matched.\n";
    }

    // determine the mode of this interaction
    //
    // Best option: compare the vertex
    if (this_vertex.match.event_vertex_id != -1) {
      this_vertex.match.mode = truth[this_vertex.match.event_vertex_id].match.mode;
    }
    // next -- try the track
    else if (track_truth != NULL) {
      if (track_truth->Origin() == simb::Origin_t::kBeamNeutrino && track_truth->GetNeutrino().CCNC() == simb::kCC) {
        this_vertex.match.mode = mCC;
      }
      else if (track_truth->Origin() == simb::Origin_t::kBeamNeutrino && track_truth->GetNeutrino().CCNC() == simb::kNC) {
        this_vertex.match.mode = mNC;
      }
      else if (track_truth->Origin() ==  simb::Origin_t::kCosmicRay) {
        this_vertex.match.mode = mCosmic;
      }
      else {
        this_vertex.match.mode = mOther;
      }
    }
    // **shrug**
    else {
      this_vertex.match.mode = mOther;
    }
    // set the completion
    this_vertex.match.track_completion = track_truth_match.completion;

    // bad vertex
    this_vertex.match.is_misreconstructed = this_vertex.match.tmode == tmNeutrino && this_vertex.match.event_vertex_id == -1;

    // get the mcturth indexes from the event indexes
    if (this_vertex.match.event_vertex_id != -1) this_vertex.match.mctruth_vertex_id = truth[this_vertex.match.event_vertex_id].match.mctruth_vertex_id;
    else this_vertex.match.mctruth_vertex_id = -1;
    if (this_vertex.match.event_track_id != -1) this_vertex.match.mctruth_track_id = truth[this_vertex.match.event_track_id].match.mctruth_track_id;
    else this_vertex.match.mctruth_track_id = -1;

    // store
    reco.push_back(std::move(this_vertex));
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

