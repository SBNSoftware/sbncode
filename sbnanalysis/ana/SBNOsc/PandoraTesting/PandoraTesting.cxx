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

#include "core/Event.hh"
#include "PandoraTesting.h"
#include "../Utilities.h"

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

#include "larcore/Geometry/Geometry.h"

namespace ana {
  namespace SBNOsc {

PandoraTesting::PandoraTesting() :
  SelectionBase(),
  _event_counter(0),
  _nu_count(0),
  _interactionInfo(new std::vector<NuMuInteraction>),
  _recoInteractionInfo(new std::vector<RecoInteractionInfo>),
  _recoParticles(new std::vector<RecoParticle>) {}

double aaBoxesMin(const std::vector<geoalgo::AABox> &boxes, unsigned dim) {
  return std::min_element(boxes.begin(), boxes.end(), [dim](auto &lhs, auto &rhs) { return lhs.Min()[dim] < rhs.Min()[dim]; })->Min()[dim];
}

double aaBoxesMax(const std::vector<geoalgo::AABox> &boxes, unsigned dim) {
  return std::max_element(boxes.begin(), boxes.end(), [dim](auto &lhs, auto &rhs) { return lhs.Max()[dim] < rhs.Max()[dim]; })->Max()[dim];
}

void PandoraTesting::Initialize(fhicl::ParameterSet* config) {
  if (config) {
    fhicl::ParameterSet pconfig = config->get<fhicl::ParameterSet>("PandoraTesting");

    // setup active volume bounding boxes
    std::vector<fhicl::ParameterSet> AVs = \
      pconfig.get<std::vector<fhicl::ParameterSet> >("active_volumes");
    for (auto const& AV : AVs) {
      double xmin = AV.get<double>("xmin");
      double ymin = AV.get<double>("ymin");
      double zmin = AV.get<double>("zmin");
      double xmax = AV.get<double>("xmax");
      double ymax = AV.get<double>("ymax");
      double zmax = AV.get<double>("zmax");
      _config.active_volumes.emplace_back(xmin, ymin, zmin, xmax, ymax, zmax);
    }

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

    _config.doFVCut = pconfig.get<bool>("doFVcut", true);
    _config.trajPointLength = pconfig.get<bool>("trajPointLength", true);
    _config.vertexDistanceCut = pconfig.get<double>("vertexDistance", -1);
    _config.minLengthContainedTrack = pconfig.get<double>("minLengthContainedTrack", -1);
    _config.minLengthExitingTrack = pconfig.get<double>("minLengthExitingTrack", -1);
    _config.trackVisibleEnergyThreshold = pconfig.get<double>("trackVisibleEnergyThreshold", 0.);
    _config.showerEnergyDistortion = pconfig.get<double>("showerEnergyDistortion", 0.);
    _config.trackEnergyDistortion = pconfig.get<double>("trackEnergyDistortion", 0.);
    _config.leptonEnergyDistortionContained = pconfig.get<double>("leptonEnergyDistortionContained", 0.);
    _config.leptonEnergyDistortionLeavingA = pconfig.get<double>("leptonEnergyDistortionLeavingA", 0.);
    _config.leptonEnergyDistortionLeavingB = pconfig.get<double>("leptonEnergyDistortionLeavingB", 0.);
    _config.acceptShakyTracks = pconfig.get<bool>("acceptShakyTracks", false);
    _config.verbose = pconfig.get<bool>("verbose", false);
    _config.cutKMEC = pconfig.get<bool>("cutKMEC", false);
    _config.onlyKMEC = pconfig.get<bool>("onlyKMEC", false);

    // setup weight config
    _config.selectionEfficiency = pconfig.get<double>("selectionEfficiency", 1.0);
    _config.uniformWeights = pconfig.get<std::vector<std::string>>("uniformWeights", {});
    _config.constantWeight = pconfig.get<double>("constantWeight", 1.0);

    _config.reco_vertex_cluster_distance = pconfig.get<double>("vertexClusterDistance", 5.0);

    // crt stuff
    _config.crt_time_limit = pconfig.get<double>("CRT_time_limit", 0.2); // from SBND CRT Track Producer
    _config.crt_avg_hit_dist = pconfig.get<double>("CRT_avg_hit_dist", 30.); // from SBND CRT Track Producer
    _config.crt_use_top_plane = pconfig.get<bool>("CRT_use_top_plane", true); // from SBND CRT Track Producer
    _config.crt_dist_limit = pconfig.get<double>("CRT_dist_limit", 25.); // from SBND CRT Track Producer

    // get tag names
    _config.HitTag = config->get<std::string>("HitTag", "gaushit");
    _config.RecoTrackTag = config->get<std::string>("RecoTrackTag", "pandoraTrack");
    _config.RecoVertexTag = config->get<std::string>("RecoVertexTag", "pandora");
    _config.PFParticleTag = config->get<std::string>("PFParticleTag", "pandora");


  }

  // Setup histo's for root output
  fOutputFile->cd();
  auto cut_names = cutNames();
  for (unsigned i = 0; i < nCuts; i++) {
    _root_histos[i].h_numu_reco_dist = new TH1D(("numu_reco_dist_" + cut_names[i]).c_str(), "numu_reco_dist", 101, -10, 1000);
    _root_histos[i].h_numu_ccqe = new TH1D(("numu_ccqe_" + cut_names[i]).c_str(), "numu_ccqe", 100, 0, 10);
    _root_histos[i].h_numu_trueE = new TH1D(("numu_trueE_" + cut_names[i]).c_str(), "numu_trueE", 100, 0 , 10);
    _root_histos[i].h_numu_visibleE = new TH1D(("numu_visibleE_" + cut_names[i]).c_str(), "numu_visibleE", 100, 0, 10);
    _root_histos[i].h_numu_true_v_visibleE = new TH1D(("numu_true_v_visibleE_" + cut_names[i]).c_str(), "numu_true_v_visibleE", 100, -10, 10);
    _root_histos[i].h_numu_t_length = new TH1D(("numu_t_length_" + cut_names[i]).c_str(), "numu_t_length", 101, -10, 1000);
    _root_histos[i].h_numu_contained_L = new TH1D(("numu_contained_L_" + cut_names[i]).c_str(), "numu_contained_L", 101, -10 , 1000);
    _root_histos[i].h_numu_t_is_contained = new TH1D(("t_is_contained_" + cut_names[i]).c_str(), "t_is_contained", 3, -1.5, 1.5);
    _root_histos[i].h_numu_t_is_muon = new TH1D(("t_is_muon_" + cut_names[i]).c_str(), "t_is_muon", 3, -1.5, 1.5);
    _root_histos[i].h_numu_Vxy = new TH2D(("numu_Vxy_" + cut_names[i]).c_str(), "numu_Vxy",
      20, aaBoxesMin(_config.active_volumes, 0), aaBoxesMax(_config.active_volumes, 0),
      20, aaBoxesMin(_config.active_volumes, 1), aaBoxesMax(_config.active_volumes, 1));
    _root_histos[i].h_numu_Vxz = new TH2D(("numu_Vxz_" + cut_names[i]).c_str(), "numu_Vxz",
      20, aaBoxesMin(_config.active_volumes, 0), aaBoxesMax(_config.active_volumes, 0),
      20, aaBoxesMin(_config.active_volumes, 2), aaBoxesMax(_config.active_volumes, 2));
    _root_histos[i].h_numu_Vyz = new TH2D(("numu_Vyz_" + cut_names[i]).c_str(), "numu_Vyz",
      20, aaBoxesMin(_config.active_volumes, 1), aaBoxesMax(_config.active_volumes, 1),
      20, aaBoxesMin(_config.active_volumes, 2), aaBoxesMax(_config.active_volumes, 2));
    _root_histos[i].h_numu_contained_L_sig = new TH1D(("numu_contained_L_sig_" + cut_names[i]).c_str(), "numu_contained_L_sig", 101, -10 , 1000);
    _root_histos[i].h_numu_contained_L_bkg = new TH1D(("numu_contained_L_bkg_" + cut_names[i]).c_str(), "numu_contained_L_bkg", 101, -10 , 1000);
    _root_histos[i].h_numu_open_angle_sig = new TH1D(("numu_open_angle_sig_" + cut_names[i]).c_str(), "numu_open_angle_sig", 60, 0., 2 * TMath::Pi());
    _root_histos[i].h_numu_open_angle_bkg = new TH1D(("numu_open_angle_bkg_" + cut_names[i]).c_str(), "numu_open_angle_bkg", 60, 0.,2 * TMath::Pi());
    _root_histos[i].h_numu_cross_TPC_sig = new TH1D(("numu_cross_TPC_sig_" + cut_names[i]).c_str(), "numu_cross_TPC_sig", 3, -1.5, 1.5);
    _root_histos[i].h_numu_cross_TPC_bkg = new TH1D(("numu_cross_TPC_bkg_" + cut_names[i]).c_str(), "numu_cross_TPC_bkg", 3, -1.5, 1.5);
    _root_histos[i].h_numu_reco_energy_sig = new TH1D(("numu_reco_energy_sig_" + cut_names[i]).c_str(), "numu_reco_energy", 50, 0, 5);
    _root_histos[i].h_numu_reco_energy_bkg = new TH1D(("numu_reco_energy_bkg_" + cut_names[i]).c_str(), "numu_reco_energy", 50, 0, 5);
  }

  // set up TGraph keeping track of cut counts
  _cut_counts = new TGraph(PandoraTesting::nCuts + 1);

  // add branches
  fTree->Branch("numu_interaction", &_interactionInfo);
  fTree->Branch("reco_interaction", &_recoInteractionInfo);
  fTree->Branch("reco_event", &_recoEventInfo);
  fTree->Branch("reco_particles", &_recoParticles);

  // setup services manager
  _manager = new core::ServiceManager(core::Detector::kSBND);

  hello();
}


void PandoraTesting::Finalize() {
  // write out histos
  fOutputFile->cd();
  for (unsigned i = 0; i < nCuts; i++) {
    _root_histos[i].h_numu_reco_dist->Write();
    _root_histos[i].h_numu_ccqe->Write();
    _root_histos[i].h_numu_trueE->Write();
    _root_histos[i].h_numu_visibleE->Write();
    _root_histos[i].h_numu_true_v_visibleE->Write();
    _root_histos[i].h_numu_t_length->Write();
    _root_histos[i].h_numu_t_is_muon->Write();
    _root_histos[i].h_numu_contained_L->Write();
    _root_histos[i].h_numu_t_is_contained->Write();
    _root_histos[i].h_numu_Vxy->Write();
    _root_histos[i].h_numu_Vxz->Write();
    _root_histos[i].h_numu_Vyz->Write();
    _root_histos[i].h_numu_contained_L_sig->Write();
    _root_histos[i].h_numu_contained_L_bkg->Write();
    _root_histos[i].h_numu_open_angle_sig->Write();
    _root_histos[i].h_numu_open_angle_bkg->Write();
    _root_histos[i].h_numu_cross_TPC_sig->Write();
    _root_histos[i].h_numu_cross_TPC_bkg->Write();
    _root_histos[i].h_numu_reco_energy_sig->Write();
    _root_histos[i].h_numu_reco_energy_bkg->Write();
  }
  _cut_counts->Write();
}


bool PandoraTesting::ProcessEvent(const gallery::Event& ev, const std::vector<Event::Interaction> &truth, std::vector<Event::RecoInteraction>& reco) {
  if (_event_counter % 10 == 0) {
    std::cout << "PandoraTesting: Processing event " << _event_counter << " "
              << "(" << _nu_count << " neutrinos selected)"
              << std::endl;
  }

  // clean up containers
  _event_counter++;
  _interactionInfo->clear();
  _recoInteractionInfo->clear();
  _recoParticles->clear();

  // Get truth
  auto const& mctruths = \
    *ev.getValidHandle<std::vector<simb::MCTruth> >(fTruthTag);
  // get tracks and showers
  auto const& mctracks = \
    *ev.getValidHandle<std::vector<sim::MCTrack> >(fMCTrackTag);
  auto const& mcshowers = \
    *ev.getValidHandle<std::vector<sim::MCShower> >(fMCShowerTag);

  // update total count of interactions
  _cut_counts->SetPoint(0, 0, _cut_counts->GetY()[0] + mctruths.size());

  // cluster the vertices in the reco info
  std::vector<RecoParticle> reco_particles = ReconstructionInfo(ev);
  _recoParticles->insert(_recoParticles->begin(), reco_particles.begin(), reco_particles.end());

  // select good vertices
  std::vector<ClusteredVertex> clustered = SelectVertices(reco_particles);

  // get reconstruction infor for the whole event
  _recoEventInfo = EventReconstructionInfo(ev, clustered);

  // Iterate through the neutrinos
  bool selected = false;
  for (size_t i=0; i<mctruths.size(); i++) {
    auto const& mctruth = mctruths.at(i);
    // get the neutrino
    const simb::MCNeutrino& nu = mctruth.GetNeutrino();

    // setup energy calculations
    VisibleEnergyCalculator calculator;
    calculator.lepton_pdgid = 13;
    calculator.track_threshold =  _config.trackVisibleEnergyThreshold;
    calculator.shower_energy_distortion = _config.showerEnergyDistortion;
    calculator.track_energy_distortion = _config.trackEnergyDistortion;
    calculator.lepton_energy_distortion_contained = _config.leptonEnergyDistortionContained;
    calculator.lepton_energy_distortion_leaving_A = _config.leptonEnergyDistortionLeavingA;
    calculator.lepton_energy_distortion_leaving_B = _config.leptonEnergyDistortionLeavingB;

    // build the interaction
    Event::Interaction interaction = truth[i];

    // Get selection-specific info
    //
    // Get the neutrino track index
    int track_ind = GetNeutrinoTrack(ev, mctruth);
    // Start with the interaction stuff
    // This also sets the lepton variables in the calculator
    NuMuInteraction intInfo = interactionInfo(ev, mctruth, track_ind, calculator);

    double visible_energy = visibleEnergy(mctruth, mctracks, mcshowers, calculator, false);

    Event::RecoInteraction reco_interaction(interaction, i);
    reco_interaction.reco_energy = visible_energy;

    // Build the weight of this event
    double weight = 1.;
    // whether this event is signal or background
    bool is_signal = abs(intInfo.t_pdgid) == 13; // muon
    // selection efficiency
    if (is_signal) {
      weight *= _config.selectionEfficiency;
    }
    //for(std::map<std::string,std::vector<double>>::iterator it = interaction.weights.begin(); it != interaction.weights.end(); ++it) {
    //  std::cout << it->first <<std::endl;
    //}
    // apply uniofrm weights (e.g. bnbcorrection)
    for (auto const &key: _config.uniformWeights) {
       weight *= interaction.weights.at(key)[0];
    }
    // apply constant weight
    weight *= _config.constantWeight;
    reco_interaction.weight = weight;

    // get the reconstruction information
    RecoInteractionInfo reco_interaction_info = Reconstruct(ev, mctruth, i, clustered, track_ind);

    // run selection
    std::array<bool, PandoraTesting::nCuts> selection = Select(ev, mctruth, i, intInfo, clustered);

    // pass iff pass each cut
    bool pass_selection = std::find(selection.begin(), selection.end(), false) == selection.end();
    if (pass_selection) {
      // store interaction in reco tree
      reco.push_back(reco_interaction);
      // store local info
      _interactionInfo->push_back(intInfo);
      _recoInteractionInfo->push_back(reco_interaction_info);

      _nu_count++;
      selected = true;
    }

    // fill histos
    for (size_t select_i=0; select_i < selection.size(); select_i++) {
      if (selection[select_i]) {
        _root_histos[select_i].h_numu_reco_dist->Fill(reco_interaction_info.reco_vertex_distance);
        _root_histos[select_i].h_numu_trueE->Fill(interaction.neutrino.energy);
        _root_histos[select_i].h_numu_ccqe->Fill(ECCQE(interaction.lepton.momentum, interaction.lepton.energy));
        _root_histos[select_i].h_numu_visibleE->Fill(visible_energy);
        _root_histos[select_i].h_numu_true_v_visibleE->Fill(visible_energy - interaction.neutrino.energy);
        _root_histos[select_i].h_numu_t_length->Fill(intInfo.t_length);
        _root_histos[select_i].h_numu_contained_L->Fill(intInfo.t_contained_length);
        _root_histos[select_i].h_numu_t_is_muon->Fill(abs(intInfo.t_pdgid) == 13);
        _root_histos[select_i].h_numu_t_is_contained->Fill(intInfo.t_is_contained);
        _root_histos[select_i].h_numu_Vxy->Fill(nu.Nu().Vx(), nu.Nu().Vy());
        _root_histos[select_i].h_numu_Vxz->Fill(nu.Nu().Vx(), nu.Nu().Vz());
        _root_histos[select_i].h_numu_Vyz->Fill(nu.Nu().Vy(), nu.Nu().Vz());
        if (isSignal(interaction)) {
          _root_histos[select_i].h_numu_contained_L_sig->Fill(intInfo.t_contained_length);
          _root_histos[select_i].h_numu_open_angle_sig->Fill(intInfo.t_open_angle);
          _root_histos[select_i].h_numu_cross_TPC_sig->Fill(intInfo.t_cross_tpc);
          _root_histos[select_i].h_numu_reco_energy_sig->Fill(visible_energy);
        }
        else { // is bkg
          _root_histos[select_i].h_numu_contained_L_bkg->Fill(intInfo.t_contained_length);
          _root_histos[select_i].h_numu_open_angle_bkg->Fill(intInfo.t_open_angle);
          _root_histos[select_i].h_numu_cross_TPC_bkg->Fill(intInfo.t_cross_tpc);
          _root_histos[select_i].h_numu_reco_energy_bkg->Fill(visible_energy);
        }

        // also update cut count
        _cut_counts->SetPoint(select_i+1, select_i+1, _cut_counts->GetY()[select_i+1] + 1);
      }
    }
  }

  return selected;
}

// get information associated with track
PandoraTesting::TrackInfo PandoraTesting::trackInfo(const sim::MCTrack &track) {
  double contained_length = 0;
  double length = 0;
  bool crosses_tpc = false;

  // If the interaction is outside the active volume, then g4 won't generate positions for the track.
  // So size == 0 => outside FV
  //
  // If size != 0, then we have to check active volume
  bool contained_in_AV = track.size() > 0;

  // Get the length and determine if any point leaves the active volume
  //
  // Use every trajectory point if configured and if the MCTrack has trajectory points
  if (track.size() != 0 && _config.trajPointLength) {
    TLorentzVector pos = track.Start().Position();
    // get the active volume that the start position is in
    int active_volume_index = -1;
    // contruct pos Point
    geoalgo::Point_t pos_point(pos);
    for (int i = 0; i < _config.active_volumes.size(); i++) {
      if (_config.active_volumes[i].Contain(pos_point)) {
        active_volume_index = i;
      }
    }

    // only consider contained length in the active volume containing the interaction
    std::vector<geoalgo::AABox> volumes;
    if (active_volume_index >= 0) {
      volumes.push_back(_config.active_volumes[active_volume_index]);
    }

    // get initial TPC of track
    geo::TPCID last_tpc_id = GetTPCIndex(pos.Vect());
    
    for (int i = 1; i < track.size(); i++) {
      // update if track is contained
      if (contained_in_AV) contained_in_AV = containedInAV(pos.Vect());
      
      // update length
      contained_length += containedLength(track[i].Position().Vect(), pos.Vect(), volumes);
      length += (track[i].Position().Vect() - pos.Vect()).Mag();

      // check if track has crossed TPC
      if (containedInAV(pos.Vect())) {
        geo::TPCID this_tpc_id = GetTPCIndex(pos.Vect());
        if (this_tpc_id && last_tpc_id && this_tpc_id != last_tpc_id) crosses_tpc = true;
        last_tpc_id = this_tpc_id;
      }
      
      pos = track[i].Position();
    }
  }
  // If active volume is misconfigured, then tracks may be generated w/out points.
  // Optionally, we can accept them.
  //
  // Also, use this method if configured
  else if (_config.acceptShakyTracks || !_config.trajPointLength) {
    contained_length = containedLength(track.Start().Position().Vect(), track.End().Position().Vect(), _config.active_volumes);
    length = (track.Start().Position().Vect() - track.End().Position().Vect()).Mag();
    contained_in_AV = containedInAV(track.Start().Position().Vect()) && containedInAV(track.End().Position().Vect());

    geo::TPCID start_tpc_id = GetTPCIndex(track.Start().Position().Vect());
    geo::TPCID end_tpc_id = GetTPCIndex(track.End().Position().Vect());
    crosses_tpc = start_tpc_id && end_tpc_id && start_tpc_id != end_tpc_id;

    //std::cout << "WARNING: SHAKY TRACK." << std::endl;
    //std::cout << "CONTAINED IN FV: " << contained_in_AV << std::endl;
    //std::cout << "CONTAINED LENGTH: " << contained_length << std::endl; 
    //std::cout << "LENGTH: " << length << std::endl;
    //
    //TVector3 start = track.Start().Position().Vect();
    //TVector3 end = track.End().Position().Vect();
    //std::cout << "START: " << start.X() << " " << start.Y() << " " << start.Z() << std::endl;
    //std::cout << "END: " << end.X() << " " << end.Y() << " " << end.Z() << std::endl;
  }
  return PandoraTesting::TrackInfo({contained_in_AV, crosses_tpc, contained_length, length});
}

int PandoraTesting::GetNeutrinoTrack(const gallery::Event &ev, const simb::MCTruth &mctruth) {
  // get handle to tracks and showers
  auto const& mctrack_list = \
    *ev.getValidHandle<std::vector<sim::MCTrack> >(fMCTrackTag);
  auto const& mcshower_list = \
    *ev.getValidHandle<std::vector<sim::MCShower> >(fMCShowerTag);

  // and particles
  auto const& mcparticle_list = \
    *ev.getValidHandle<std::vector<simb::MCParticle>>(fMCParticleTag);

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

PandoraTesting::NuMuInteraction PandoraTesting::interactionInfo(const gallery::Event &ev, const simb::MCTruth &mctruth, int track_ind, VisibleEnergyCalculator &calculator) {
  auto const& mctrack_list = \
    *ev.getValidHandle<std::vector<sim::MCTrack> >(fMCTrackTag);

  // if there's no track, return nonsense
  if (track_ind == -1) {
    // set calculator variables
    calculator.lepton_contained = false;
    calculator.lepton_contained_length = -1;
    calculator.lepton_index = -1;
    return PandoraTesting::NuMuInteraction({false, false, -1, -1, -1, -1, -1, -1});
  }
  // otherwise get the track info and energy info
  else {
    // get track info for our lepton
    PandoraTesting::TrackInfo t_info = trackInfo(mctrack_list[track_ind]);

    // set calculator variables
    calculator.lepton_contained = t_info.t_is_contained; 
    calculator.lepton_contained_length = t_info.t_contained_length;
    calculator.lepton_index = track_ind;

    // smear the energy
    double smeared_energy = smearLeptonEnergy(mctrack_list[track_ind], calculator);
    // truth kinetic energy
    double truth_energy = (mctrack_list[track_ind].Start().E()) / 1000.; /* MeV -> GeV */
    // opening angle
    double open_angle = mctruth.GetNeutrino().Nu().Momentum().Vect().Angle(mctrack_list[track_ind].Start().Momentum().Vect());

    return PandoraTesting::NuMuInteraction({t_info.t_is_contained, t_info.t_cross_tpc, t_info.t_contained_length, t_info.t_length, mctrack_list[track_ind].PdgCode(), truth_energy, smeared_energy, open_angle});
  }
}

std::vector<PandoraTesting::ClusteredVertex> PandoraTesting::SelectVertices(const std::vector<PandoraTesting::RecoParticle>& reco_particles) {
  std::vector<PandoraTesting::ClusteredVertex> ret;
  for (auto const &particle: reco_particles) {
     if (particle.p_is_neutrino) {
       ret.insert(ret.begin(), particle.vertices.begin(), particle.vertices.end());
     }
  }
  return ret;
}

std::vector<PandoraTesting::RecoParticle> PandoraTesting::ReconstructionInfo(const gallery::Event &event) {
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
  std::vector<PandoraTesting::RecoParticle> ret;

  // iterate over all of the pfparticles
  for (size_t i = 0; i < pfp_handle->size(); i++) {
    // new reco particle
    PandoraTesting::RecoParticle this_particle;

    // get the PFParticle
    const recob::PFParticle& this_pfp = pfp_handle->at(i);
    // get the metadata
    const larpandoraobj::PFParticleMetadata* this_metadata = pfp_metadatas.at(i).at(0);
    // and the properties dict
    auto const &properties = this_metadata->GetPropertiesMap();
    // get the reco vertices
    const std::vector<const recob::Vertex*> &this_vertices = pfp_vertices.at(i);

    // cluster the vertices
    this_particle.vertices = ClusterVertices(this_vertices);

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

std::vector<PandoraTesting::ClusteredVertex> PandoraTesting::ClusterVertices(const std::vector<const recob::Vertex*> &reco_vertices) {

  // returned
  std::vector<ClusteredVertex> ret;

  // cluster all vertices within a given length

  // already clustered vertices
  std::set<size_t> clustered;

  // TODO: make more efficient?
  //
  // run until all are clustered
  while (clustered.size() < reco_vertices.size()) {
    for (size_t vertex_i = 0; vertex_i < reco_vertices.size(); vertex_i++) {
      if (!clustered.count(vertex_i)) {
        clustered.insert(vertex_i);
        // found a new vertex to cluster! Make a new cluster
        ClusteredVertex vert;
        // setup this vertex
        auto pos = reco_vertices[vertex_i]->position();
        vert.vertices.emplace_back(pos.X(), pos.Y(), pos.Z());

        // seed with a single vertex, re-run until no more are pulled in
        bool added = true;
        while (added) {
          added = false;
          for (size_t vertex_j = 0; vertex_j < reco_vertices.size(); vertex_j++) {
            if (!clustered.count(vertex_j)) {
              auto this_pos = reco_vertices[vertex_j]->position();
              TVector3 this_vertex(this_pos.X(), this_pos.Y(), this_pos.Z());
              bool do_cluster = false;
              for (TVector3 &vertex: vert.vertices) {
                if ((this_vertex - vertex).Mag() < _config.reco_vertex_cluster_distance) {
                  do_cluster = true;
                  break;
                }
              }
              if (do_cluster) {
                clustered.insert(vertex_j);
                vert.vertices.push_back(this_vertex);
                added = true;
              }
            }
          }
        }
        vert.CalcMean();
        ret.push_back(vert);
      }
    }
  }
  return ret;
}

PandoraTesting::RecoEventInfo PandoraTesting::EventReconstructionInfo(const gallery::Event &ev, std::vector<ClusteredVertex> &vertices) {
  RecoEventInfo ret {0, 0};
  for (const ClusteredVertex &v: vertices) {
    ret.n_vertices += v.vertices.size();
    ret.n_clustered_vertices += 1;
  }
  return ret;
}

PandoraTesting::RecoInteractionInfo PandoraTesting::Reconstruct(const gallery::Event &ev, const simb::MCTruth &mctruth, unsigned truth_ind, std::vector<PandoraTesting::ClusteredVertex> &vertices, int track_ind) {
  // get the hits
  auto const& gaus_hits = \
    *ev.getValidHandle<std::vector<recob::Hit>>(_config.HitTag);

  // get the reco tracks
  auto const &reco_tracks = \
    *ev.getValidHandle<std::vector<recob::Track>>(_config.RecoTrackTag);

  // and the truth tracks
  auto const& mctrack_list = \
    *ev.getValidHandle<std::vector<sim::MCTrack> >(fMCTrackTag);

  // get the location of this neutrino
  TVector3 truth_vertex =  mctruth.GetNeutrino().Nu().Trajectory().Position(0).Vect();

  // get the distance to the closest vertex
  double distance_to_reco_vertex = -1;
  for (const ClusteredVertex &v: vertices) {
    double this_distance = (truth_vertex - v.mean).Mag();
    if (distance_to_reco_vertex < 0 || distance_to_reco_vertex > this_distance) {
      distance_to_reco_vertex = this_distance;
    }
  }

  // get the number of vertices within tolerance if there is a valid track
  int n_track_clustered_vertices = 0;
  if (track_ind >= 0 && mctrack_list[track_ind].size() > 0) {
    for (const ClusteredVertex &v: vertices) {
      const sim::MCTrack &track = mctrack_list[track_ind];
      TVector3 pos = track[0].Position().Vect();
      for (int i = 1; i < track.size(); i++) {
        TVector3 this_pos = track[i].Position().Vect();
        if (closestDistance(pos, this_pos, v.mean) < 5.0) {
          n_track_clustered_vertices += 1;
          break;
        }
      }
    }
  }
  else {
    n_track_clustered_vertices = -1;
  }

  RecoInteractionInfo ret;
  ret.reco_vertex_distance = distance_to_reco_vertex;
  ret.n_track_clustered_vertices = n_track_clustered_vertices;

  return ret;
}

std::array<bool, PandoraTesting::nCuts> PandoraTesting::Select(const gallery::Event& ev, const simb::MCTruth& mctruth,
      unsigned truth_ind, const PandoraTesting::NuMuInteraction &intInfo, const std::vector<ClusteredVertex> &reco_vertices) {
  // get the neutrino
  const simb::MCNeutrino& nu = mctruth.GetNeutrino();

  // has valid track
  bool pass_valid_track = intInfo.t_pdgid != -1;

  // pass fiducial volume cut
  bool pass_FV = passFV(nu.Nu().Position().Vect());

  // min length cut
  bool pass_min_length = passMinLength(intInfo.t_contained_length, intInfo.t_is_contained);

  // pass vertex reconstruction cut
  bool pass_reco_vertex = passRecoVertex(nu.Nu().Position().Vect(), reco_vertices);

  // print selection information
  if (_config.verbose) {
    std::cout << "NEW EVENT" << std::endl;
    std::cout << "CCNC: " << nu.CCNC() << " MODE: " << nu.Mode() << " PDG: " << nu.Nu().PdgCode() << std::endl;
    std::cout << "Track PDG: " << intInfo.t_pdgid <<std::endl;
    std::cout << "pass Valid Track: " << pass_valid_track << std::endl;
    std::cout << "Pos: " << nu.Nu().Vx() << " " << nu.Nu().Vy() << " " << nu.Nu().Vz() << std::endl;
    std::cout << "pass FV: " << pass_FV << std::endl;
    std::cout << "pass Reco: " << pass_reco_vertex << std::endl;
    std::cout << "Length: " << intInfo.t_contained_length << " Contained: " << intInfo.t_is_contained << std::endl;
    std::cout << "pass Length: " << pass_min_length << std::endl;

    std::cout << "Vertex: " << nu.Nu().Position().Vect().X() << " " << nu.Nu().Position().Vect().Y() << " " << nu.Nu().Position().Vect().Z() << std::endl;
    for (auto const &v: reco_vertices) {
      std::cout << "Clustered vertex: " << v.mean.X() << " " << v.mean.Y() << " " << v.mean.Z() << std::endl;
      for (auto const &vv: v.vertices) {
        std::cout << "Component vertex: " << vv.X() << " " << vv.Y() << " " << vv.Z() << std::endl;
      }
    } 
  }

  // TEMP: add in an active volume cut
  bool pass_AV = containedInAV(nu.Nu().Position().Vect()); 

  // STUDY KMEC: remove MEC events
  bool pass_kMEC = !(_config.cutKMEC && nu.Mode() == simb::kMEC) && !(_config.onlyKMEC && nu.Mode() != simb::kMEC);

  // retrun list of cuts
  return {
    pass_kMEC, 
    pass_AV && pass_kMEC,
    pass_valid_track && pass_kMEC && pass_AV,
    pass_valid_track && pass_kMEC && pass_AV && pass_FV,
    pass_valid_track && pass_kMEC && pass_AV && pass_FV && pass_min_length,
    pass_valid_track && pass_kMEC && pass_AV && pass_FV && pass_min_length && pass_reco_vertex
  };
}

bool PandoraTesting::containedInAV(const TVector3 &v) {
  geoalgo::Point_t p(v); 
  for (auto const& AV: _config.active_volumes) {
    if (AV.Contain(p)) return true;
  }
  return false;
}

bool PandoraTesting::containedInFV(const TVector3 &v) {
  geoalgo::Point_t p(v);
  for (auto const& FV: _config.fiducial_volumes) {
    if (FV.Contain(p)) return true;
  }
  return false;
}

bool PandoraTesting::passRecoVertex(const TVector3 &truth_v, const std::vector<ClusteredVertex> &vertices) {
  if (_config.vertexDistanceCut < 0) return true;

  for (auto const &v: vertices) {
    if ((truth_v - v.mean).Mag() < _config.vertexDistanceCut) {
      return true;
    }
  }

  return false;
}

bool PandoraTesting::passMinLength(double length, bool stop_in_tpc) {
  if (!stop_in_tpc)
    return _config.minLengthExitingTrack < 0 || length > _config.minLengthExitingTrack;
  else
    return _config.minLengthContainedTrack < 0 || length > _config.minLengthContainedTrack;
}

  }  // namespace SBNOsc
}  // namespace ana


DECLARE_SBN_PROCESSOR(ana::SBNOsc::PandoraTesting)

