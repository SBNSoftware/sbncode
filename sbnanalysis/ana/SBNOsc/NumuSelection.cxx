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

#include "core/Event.hh"
#include "NumuSelection.h"
#include "Utilities.h"

namespace ana {
  namespace SBNOsc {

NumuSelection::NumuSelection() :
  SelectionBase(),
  _event_counter(0),
  _nu_count(0),
  _interactionInfo(new std::vector<NuMuInteraction>) {}

double aaBoxesMin(const std::vector<geoalgo::AABox> &boxes, unsigned dim) {
  return std::min_element(boxes.begin(), boxes.end(), [dim](auto &lhs, auto &rhs) { return lhs.Min()[dim] < rhs.Min()[dim]; })->Min()[dim];
}

double aaBoxesMax(const std::vector<geoalgo::AABox> &boxes, unsigned dim) {
  return std::max_element(boxes.begin(), boxes.end(), [dim](auto &lhs, auto &rhs) { return lhs.Max()[dim] < rhs.Max()[dim]; })->Max()[dim];
}

void NumuSelection::Initialize(fhicl::ParameterSet* config) {
  if (config) {
    fhicl::ParameterSet pconfig = config->get<fhicl::ParameterSet>("NumuSelection");

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
    _config.selectMode = pconfig.get<int>("selectMode", -1);
    _config.selectCCNC = pconfig.get<int>("selectCCNC", -1);
    _config.onlyKMEC = pconfig.get<bool>("onlyKMEC", false);

    // setup weight config
    _config.selectionEfficiency = pconfig.get<double>("selectionEfficiency", 1.0);
    _config.uniformWeights = pconfig.get<std::vector<std::string>>("uniformWeights", {});
    _config.constantWeight = pconfig.get<double>("constantWeight", 1.0);

  }

  // Setup histo's for root output
  fOutputFile->cd();
  auto cut_names = cutNames();
  for (unsigned i = 0; i < nCuts; i++) {
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

    _root_histos[i].h_numu_Vx_sig = new TH1D(("numu_Vx_sig" + cut_names[i]).c_str(), "numu_Vx_sig", 20, 
      aaBoxesMin(_config.fiducial_volumes, 0), aaBoxesMax(_config.fiducial_volumes, 0));
    _root_histos[i].h_numu_Vy_sig = new TH1D(("numu_Vy_sig" + cut_names[i]).c_str(), "numu_Vy_sig", 20, 
      aaBoxesMin(_config.fiducial_volumes, 1), aaBoxesMax(_config.fiducial_volumes, 1));
    _root_histos[i].h_numu_Vz_sig = new TH1D(("numu_Vz_sig" + cut_names[i]).c_str(), "numu_Vz_sig", 20, 
      aaBoxesMin(_config.fiducial_volumes, 2), aaBoxesMax(_config.fiducial_volumes, 2));

    _root_histos[i].h_numu_Vx_bkg = new TH1D(("numu_Vx_bkg" + cut_names[i]).c_str(), "numu_Vx_bkg", 20, 
      aaBoxesMin(_config.fiducial_volumes, 0), aaBoxesMax(_config.fiducial_volumes, 0));
    _root_histos[i].h_numu_Vy_bkg = new TH1D(("numu_Vy_bkg" + cut_names[i]).c_str(), "numu_Vy_bkg", 20, 
      aaBoxesMin(_config.fiducial_volumes, 1), aaBoxesMax(_config.fiducial_volumes, 1));
    _root_histos[i].h_numu_Vz_bkg = new TH1D(("numu_Vz_bkg" + cut_names[i]).c_str(), "numu_Vz_bkg", 20, 
      aaBoxesMin(_config.fiducial_volumes, 2), aaBoxesMax(_config.fiducial_volumes, 2));

    _root_histos[i].h_numu_t_is_muon_sig = new TH1D(("t_is_muon_sig_"  + cut_names[i]).c_str(), "t_is_muon_sig", 3, -1.5, 1.5);
    _root_histos[i].h_numu_t_is_muon_bkg = new TH1D(("t_is_muon_bkg_"  + cut_names[i]).c_str(), "t_is_muon_bkg", 3, -1.5, 1.5);

  }

  // set up TGraph keeping track of cut counts
  _cut_counts = new TGraph(NumuSelection::nCuts + 1);

  // add branches
  fTree->Branch("numu_interaction", &_interactionInfo);

  hello();
}


void NumuSelection::Finalize() {
  // write out histos
  fOutputFile->cd();
  for (unsigned i = 0; i < nCuts; i++) {
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

    _root_histos[i].h_numu_Vx_sig->Write();
    _root_histos[i].h_numu_Vy_sig->Write();
    _root_histos[i].h_numu_Vz_sig->Write();
    _root_histos[i].h_numu_Vx_bkg->Write();
    _root_histos[i].h_numu_Vy_bkg->Write();
    _root_histos[i].h_numu_Vz_bkg->Write();

    _root_histos[i].h_numu_t_is_muon_sig->Write();
    _root_histos[i].h_numu_t_is_muon_bkg->Write();
  }
  _cut_counts->Write();
}


bool NumuSelection::ProcessEvent(const gallery::Event& ev, const std::vector<Event::Interaction> &truth, std::vector<Event::RecoInteraction>& reco) {
  if (_event_counter % 10 == 0) {
    std::cout << "NumuSelection: Processing event " << _event_counter << " "
              << "(" << _nu_count << " neutrinos selected)"
              << std::endl;
  }

  // clean up containers
  _event_counter++;
  _interactionInfo->clear();

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
    // Start with the interaction stuff
    // This also sets the lepton variables in the calculator
    NuMuInteraction intInfo = interactionInfo(ev, mctruth, calculator);

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
    // apply uniofrm weights (e.g. bnbcorrection)
    for (auto const &key: _config.uniformWeights) {
       weight *= interaction.weights.at(key)[0];
    }
    // apply constant weight
    weight *= _config.constantWeight;
    reco_interaction.weight = weight;

    // run selection
    std::array<bool, NumuSelection::nCuts> selection = Select(ev, mctruth, i, intInfo);

    // pass iff pass each cut
    bool pass_selection = std::find(selection.begin(), selection.end(), false) == selection.end();
    if (pass_selection) {
      // store interaction in reco tree
      reco.push_back(reco_interaction);
      // store local info
      _interactionInfo->push_back(intInfo);

      _nu_count++;
      selected = true;
    }

    // fill histos
    for (size_t select_i=0; select_i < selection.size(); select_i++) {
      if (selection[select_i]) {
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

        if (is_signal) {
          _root_histos[select_i].h_numu_Vx_sig->Fill(nu.Nu().Vx()); 
          _root_histos[select_i].h_numu_Vy_sig->Fill(nu.Nu().Vy()); 
          _root_histos[select_i].h_numu_Vz_sig->Fill(nu.Nu().Vz()); 

          _root_histos[select_i].h_numu_t_is_muon_sig->Fill(abs(intInfo.t_pdgid) == 13);
        }
        else {
          _root_histos[select_i].h_numu_Vx_bkg->Fill(nu.Nu().Vx()); 
          _root_histos[select_i].h_numu_Vy_bkg->Fill(nu.Nu().Vy()); 
          _root_histos[select_i].h_numu_Vz_bkg->Fill(nu.Nu().Vz()); 

          _root_histos[select_i].h_numu_t_is_muon_bkg->Fill(abs(intInfo.t_pdgid) == 13);
        }

        // also update cut count
        _cut_counts->SetPoint(select_i+1, select_i+1, _cut_counts->GetY()[select_i+1] + 1);
      }
    }
  }

  return selected;
}

// get information associated with track
NumuSelection::TrackInfo NumuSelection::trackInfo(const sim::MCTrack &track) {
  double contained_length = 0;
  double length = 0;
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
    
    for (int i = 1; i < track.size(); i++) {
      // update if track is contained
      if (contained_in_AV) contained_in_AV = containedInAV(pos.Vect());
      
      // update length
      contained_length += containedLength(track[i].Position().Vect(), pos.Vect(), volumes);
      length += (track[i].Position().Vect() - pos.Vect()).Mag();
      
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
  return NumuSelection::TrackInfo({contained_in_AV, contained_length, length});
}

NumuSelection::NuMuInteraction NumuSelection::interactionInfo(const gallery::Event &ev, const simb::MCTruth &mctruth, VisibleEnergyCalculator &calculator) {
  // get handle to tracks and showers
  auto const& mctrack_list = \
    *ev.getValidHandle<std::vector<sim::MCTrack> >(fMCTrackTag);
  auto const& mcshower_list = \
    *ev.getValidHandle<std::vector<sim::MCShower> >(fMCShowerTag);

  // print out track/shower info
  if (_config.verbose) {
    std::cout << "\n\nINTERACTION:\n";
    std::cout << "MODE: " << mctruth.GetNeutrino().Mode() << std::endl;
    std::cout << "CC: " << mctruth.GetNeutrino().CCNC() << std::endl;
    std::cout << "\nTRACKS:\n";
    for (auto const &mct: mctrack_list) {
      std::cout << "TRACK:\n";
      std::cout << "PDG: " << mct.PdgCode() << std::endl;
      std::cout << "Vertex: " << isFromNuVertex(mctruth, mct) << std::endl;
      std::cout << "Energy: " << mct.Start().E() << std::endl;
      std::cout << "Kinetic: " << mct.Start().E() - PDGMass(mct.PdgCode()) << std::endl;
      std::cout << "Length: " << (mct.Start().Position().Vect() - mct.End().Position().Vect()).Mag() << std::endl;
      std::cout << "Process: " << mct.Process() << std::endl;
    }
    std::cout << "\nSHOWERS:\n";
    for (auto const &mcs: mcshower_list) {
      std::cout << "SHOWER:\n";
      std::cout << "PDG: " << mcs.PdgCode() << std::endl;
      std::cout << "Vertex: " << isFromNuVertex(mctruth, mcs) << std::endl;
      std::cout << "Energy: " << mcs.Start().E() << std::endl;
      std::cout << "Kinetic: " << mcs.Start().E() - PDGMass(mcs.PdgCode()) << std::endl;
      std::cout << "Length: " << (mcs.Start().Position().Vect() - mcs.End().Position().Vect()).Mag() << std::endl;
      std::cout << "Process: " << mcs.Process() << std::endl;
    }
  }

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
        //double this_contained_length = trackInfo(mctrack_list[i]).t_contained_length; 
        //if (track_contained_length < 0 || this_contained_length > track_contained_length) {
        //  track_ind = i;
        // track_contained_length = this_contained_length;
        //}
      }
    }
  }

  // if there's no track, return nonsense
  if (track_ind == -1) {
    // set calculator variables
    calculator.lepton_contained = false;
    calculator.lepton_contained_length = -1;
    calculator.lepton_index = -1;
    return NumuSelection::NuMuInteraction({false, -1, -1, -1, -1, -1});
  }
  // otherwise get the track info and energy info
  else {
    // get track info for our lepton
    NumuSelection::TrackInfo t_info = trackInfo(mctrack_list[track_ind]);

    // set calculator variables
    calculator.lepton_contained = t_info.t_is_contained; 
    calculator.lepton_contained_length = t_info.t_contained_length;
    calculator.lepton_index = track_ind;

    // smear the energy
    double smeared_energy = smearLeptonEnergy(mctrack_list[track_ind], calculator);
    // truth kinetic energy
    double truth_energy = (mctrack_list[track_ind].Start().E()) / 1000.; /* MeV -> GeV */
    return NumuSelection::NuMuInteraction({t_info.t_is_contained, t_info.t_contained_length, t_info.t_length, mctrack_list[track_ind].PdgCode(), truth_energy, smeared_energy});
  }
}

std::array<bool, NumuSelection::nCuts> NumuSelection::Select(const gallery::Event& ev, const simb::MCTruth& mctruth,
      unsigned truth_ind, const NumuSelection::NuMuInteraction &intInfo) {
  // get the neutrino
  const simb::MCNeutrino& nu = mctruth.GetNeutrino();

  // has valid track
  bool pass_valid_track = intInfo.t_pdgid != -1;

  // pass fiducial volume cut
  bool pass_FV = passFV(nu.Nu().Position().Vect());

  // min length cut
  bool pass_min_length = passMinLength(intInfo.t_contained_length, intInfo.t_is_contained);

  // pass vertex reconstruction cut
  bool pass_reco_vertex = true;
  if (_config.vertexDistanceCut > 0) {
    double truth_v[3];
    truth_v[0] = nu.Nu().Vx();
    truth_v[1] = nu.Nu().Vy();
    truth_v[2] = nu.Nu().Vz();

    // get the reco vertex information
    auto const pfp_handle = ev.getValidHandle<std::vector<recob::PFParticle>>("pandoraNu");
    art::FindMany<recob::Vertex> fvtx(pfp_handle, ev, "pandoraNu");
    // index of truth info is same as index of vertex info
    std::vector<const recob::Vertex*> vertices = fvtx.at(truth_ind);
    // neutrino will only have one associated vetex
    auto vertex = vertices[0]->position();
    TVector3 reco_v(vertex.X(), vertex.Y(), vertex.Z());

    pass_reco_vertex = passRecoVertex(nu.Nu().Position().Vect(), reco_v);
  }

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
  }

  // TEMP: add in an active volume cut
  bool pass_AV = false; 
  geoalgo::Point_t interaction(nu.Nu().Position().Vect());
  for (auto const& AV: _config.active_volumes) {
    if (AV.Contain(interaction)) pass_AV = true;
  }

  // STUDY KMEC: remove MEC events
  bool pass_kMEC = !(_config.cutKMEC && nu.Mode() == simb::kMEC) && !(_config.onlyKMEC && nu.Mode() != simb::kMEC);
  // select another mode if necessary
  bool pass_Mode = _config.selectMode < 0 || nu.Mode() == _config.selectMode;
  // maybe require cc or nc
  bool pass_CCNC = _config.selectCCNC < 0 || nu.CCNC() == _config.selectCCNC;
  pass_kMEC = pass_kMEC && pass_Mode && pass_CCNC;

  // retrun list of cuts
  return {pass_kMEC, pass_AV && pass_kMEC, pass_valid_track && pass_kMEC && pass_AV, pass_valid_track && pass_kMEC && pass_FV, pass_valid_track && pass_kMEC && pass_FV && pass_min_length};
}

bool NumuSelection::containedInAV(const TVector3 &v) {
  geoalgo::Point_t p(v); 
  for (auto const& AV: _config.active_volumes) {
    if (AV.Contain(p)) return true;
  }
  return false;
}

bool NumuSelection::containedInFV(const TVector3 &v) {
  geoalgo::Point_t p(v);
  for (auto const& FV: _config.fiducial_volumes) {
    if (FV.Contain(p)) return true;
  }
  return false;
}

bool NumuSelection::passRecoVertex(const TVector3 &truth_v, const TVector3 &reco_v) {
  if (_config.vertexDistanceCut < 0) return true;

  return (truth_v - reco_v).Mag() < _config.vertexDistanceCut;
}

bool NumuSelection::passMinLength(double length, bool stop_in_tpc) {
  if (!stop_in_tpc)
    return _config.minLengthExitingTrack < 0 || length > _config.minLengthExitingTrack;
  else
    return _config.minLengthContainedTrack < 0 || length > _config.minLengthContainedTrack;
}

  }  // namespace SBNOsc
}  // namespace ana


DECLARE_SBN_PROCESSOR(ana::SBNOsc::NumuSelection)

