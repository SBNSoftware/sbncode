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
#include "nusimdata/SimulationBase/MCTruth.h"
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
#include "NumuRecoSelection.h"
#include "../SBNOsc/Utilities.h"

#include "ubcore/LLBasicTool/GeoAlgo/GeoAABox.h"

namespace ana {
  namespace SBNOsc {

NumuRecoSelection::NumuRecoSelection() :
  SelectionBase(),
  _event_counter(0),
  _nu_count(0) {}

void NumuRecoSelection::Initialize(fhicl::ParameterSet* config) {
  _manager = new core::ServiceManager(core::Detector::kSBND);
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
      _config.fiducial_volumes.emplace_back(xmin, ymin, zmin, xmax, ymax, zmax);
    }

    // get the active and cryo volumes
    for (auto const &cryo: _manager->GetGeometryService()->IterateCryostats()) {
      _config.cryostat_volumes.push_back(cryo.BoundingBox());
      std::vector<geo::BoxBoundedGeo> this_tpc_volumes;

      geo::GeometryCore::TPC_iterator iTPC = _manager->GetGeometryService()->begin_TPC(cryo.ID()),
                                      tend = _manager->GetGeometryService()->end_TPC(cryo.ID());
      while (iTPC != tend) {
        geo::TPCGeo const& TPC = *iTPC;
        this_tpc_volumes.push_back(TPC.ActiveBoundingBox());
        iTPC++;
      }
     _config.tpc_volumes.push_back(std::move(this_tpc_volumes));
    }
   

    _config.shakyMCTracks = pconfig.get<bool>("shakyMCTracks", false);
    _config.verbose = pconfig.get<bool>("verbose", false);

    _config.requireMatched = pconfig.get<bool>("requireMatched", false);
    _config.requireContained = pconfig.get<bool>("requireContained", false);

    // setup weight config
    _config.uniformWeights = pconfig.get<std::vector<std::string>>("uniformWeights", {});
    _config.constantWeight = pconfig.get<double>("constantWeight", 1.0);

    // get tag names
    _config.HitTag = config->get<std::string>("HitTag", "gaushit");
    _config.RecoTrackTag = config->get<std::string>("RecoTrackTag", "pandoraTrack");
    _config.RecoVertexTag = config->get<std::string>("RecoVertexTag", "pandora");
    _config.PFParticleTag = config->get<std::string>("PFParticleTag", "pandora");
  }

  // Setup histo's for root output
  fOutputFile->cd();
  unsigned histo_ind = 0;
  for (unsigned i = 0; i < NumuRecoSelection::nCuts; i++) {
    for (const auto mode: allModes) {
      _root_histos[histo_ind].track_length = new TH1D(("track_length_" + mode2Str(mode) + "_" + cutNames[i]).c_str(), "track_length", -10, 1000, 101);
      histo_ind ++;
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
  for (unsigned i = 0; i < NumuRecoSelection::nCuts*NumuRecoSelection::nModes; i++) {
    _root_histos[i].track_length->Write();
  }
}

Event::RecoInteraction NumuRecoSelection::CoreRecoInteraction(const std::vector<Event::Interaction> &truth, const NumuRecoSelection::RecoVertex &vertex, double weight) {
  Event::RecoInteraction ret;
  if (vertex.match > 0) {
    ret.truth_index = vertex.match;
    ret.truth = truth[vertex.match];
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
    std::array<bool, NumuRecoSelection::nRecoCuts> cuts = ProcessRecoCuts(_recoEvent, i);
    bool pass_selection = SelectReco(cuts);

    // compute the weight for this interaction
    double weight = 1.;
    weight *= _config.constantWeight;
    // TODO: what about cosmics?
    if (vertex.match >= 0) {
      for (auto const &key: _config.uniformWeights) {
         weight *= core_truth[vertex.match].weights.at(key)[0];
      }
    }
    if (vertex.mode == mCosmic) {
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
    for (size_t cut_i=0; cut_i < NumuRecoSelection::nRecoCuts; cut_i++) {
      if (cuts[cut_i]) {
	unsigned mode_index = (NumuRecoSelection::nTruthCuts + cut_i)*NumuRecoSelection::nModes + NumuRecoSelection::nTruthCuts + vertex.mode;
	unsigned mall_index = (NumuRecoSelection::nTruthCuts + cut_i)*NumuRecoSelection::nModes + NumuRecoSelection::nTruthCuts + mAll;
	_root_histos[mode_index].track_length->Fill(vertex.track.length);
	_root_histos[mall_index].track_length->Fill(vertex.track.length);
      }
    }
  }
  // fill histos for truth information
  for (unsigned i = 0; i < _recoEvent.truth.size(); i++) {
    const NumuRecoSelection::RecoVertex &vertex = _recoEvent.truth[i];

    std::array<bool, NumuRecoSelection::nTruthCuts> cuts = ProcessTruthCuts(_recoEvent, i);

    for (size_t cut_i = 0; cut_i < NumuRecoSelection::nTruthCuts; cut_i++) {
      if (cuts[cut_i]) {
	unsigned mode_index = cut_i*NumuRecoSelection::nModes + vertex.mode;
	unsigned mall_index = cut_i*NumuRecoSelection::nModes + mAll;
	_root_histos[mode_index].track_length->Fill(vertex.track.length);
	_root_histos[mall_index].track_length->Fill(vertex.track.length);
      }
    }
  }

  return selected;
}

std::array<bool, NumuRecoSelection::nRecoCuts> NumuRecoSelection::ProcessRecoCuts(const RecoEvent &event, unsigned reco_vertex_index) {
  bool is_reco = true;
  bool is_matched = event.reco[reco_vertex_index].match >= 0;
  bool is_contained = event.reco[reco_vertex_index].track.contained_in_cryo;
  return {
    is_reco,
    is_matched,
    is_contained
  };
}

std::array<bool, NumuRecoSelection::nTruthCuts> NumuRecoSelection::ProcessTruthCuts(const RecoEvent &event, unsigned truth_vertex_index) {
  bool is_truth = true;
  bool has_track = event.truth[truth_vertex_index].track.length > 0;
  bool is_matched = event.truth[truth_vertex_index].match >= 0;
  return {
    is_truth,
    has_track,
    is_matched
  };
}

bool NumuRecoSelection::SelectReco(std::array<bool, NumuRecoSelection::nRecoCuts> &cuts) {
  return 
    cuts[0] && 
    (cuts[1] || !_config.requireMatched) &&
    (cuts[2] || !_config.requireContained);
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
  double open_angle = truth.GetNeutrino().Nu().Momentum().Vect().Angle(track.Start().Momentum().Vect());
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
    open_angle,
    pdgid,
    contained_in_cryo,
    contained_in_tpc,
    crosses_tpc
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
  for (auto const &truth: mctruths) {
    auto neutrino = truth.GetNeutrino();
    TVector3 position = neutrino.Nu().Position().Vect();

    if (!containedInFV(position)) continue;

    // get the info
    RecoVertex this_interaction;
    if (neutrino.CCNC() == simb::kCC) {
      this_interaction.mode = mCC;
    }
    else {
      this_interaction.mode = mNC;
    }
    this_interaction.position = position;
    this_interaction.nu_energy = neutrino.Nu().E();
    this_interaction.track = MCTruthTrackInfo(event, truth);
    // no match by default -- Reconstruction will set this if necessary
    this_interaction.match = -1;
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


// TODO: Implement
NumuRecoSelection::TrackInfo NumuRecoSelection::RecoTrackInfo(const std::vector<NumuRecoSelection::RecoParticle> &reco_particles, const NumuRecoSelection::RecoParticle& vertex) {
  return {-1, -1, -1, -1, false, false, false};
}

std::vector<NumuRecoSelection::RecoParticle> NumuRecoSelection::SelectVertices(const std::vector<NumuRecoSelection::RecoParticle>& reco_particles) {
  std::vector<NumuRecoSelection::RecoParticle> ret;
  for (auto const &particle: reco_particles) {
     if (particle.p_is_neutrino) {
       ret.push_back(particle);
     }
  }
  return ret;
}

NumuRecoSelection::RecoEvent NumuRecoSelection::Reconstruct(const gallery::Event &ev, std::vector<NumuRecoSelection::RecoVertex> truth) {
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

    // Try to match against truth events
    int truth_match = -1;
    for (auto const &vert: truth) {
      for (int truth_i = 0; truth_i < truth.size(); truth_i++) {
        // found a match!
        if ((truth[truth_i].position - reco_position).Mag() < 5.0) {
          truth_match = truth_i;
        }
      }
    }
    // if you found a match, set the mode
    if (truth_match >= 0) {
      assert(truth[truth_match].match == -1);
      truth[truth_match].match = reco_i;

      this_vertex.match = truth_match;
      this_vertex.mode = truth[truth_match].mode;
      
    }
    else {
      this_vertex.mode = mCosmic;
      this_vertex.match = -1;
    }

    // get the reco track info
    this_vertex.track = RecoTrackInfo(reco_particles, particle);
    // TODO: get the enrgy
    this_vertex.nu_energy = -1;
  
    // store
    reco.push_back(std::move(this_vertex));
  }

  RecoEvent event;
  event.truth = std::move(truth);
  event.reco = std::move(reco);
  return std::move(event);
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

