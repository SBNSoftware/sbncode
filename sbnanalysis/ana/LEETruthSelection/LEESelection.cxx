#include <cassert>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "uboone/EventWeight/MCEventWeight.h"
#include "gallery/Event.h"

#include <TLorentzVector.h>
#include <TRandom.h>

#include "core/SelectionBase.hh"
#include "Util.h"
#include "Config.h"
#include "Cuts.h"
#include "LEESelection.h"

namespace ana {
  namespace LEETruthSelection {

void LEESelection::Initialize(Json::Value* config) {
  // Load configuration
  fConfig.Initialize(config);

  // Add branches to the output tree
  AddBranch("np", &fOutputData.np);
  AddBranch("ntrk", &fOutputData.ntrk);
  AddBranch("dataset", &fOutputData.dataset);
  AddBranch("bnbweight", &fOutputData.bnbweight);

  AddBranch("lpid", &fOutputData.lpid);
  AddBranch("lpdg", &fOutputData.lpdg);
  AddBranch("lexit", &fOutputData.lexit);
  AddBranch("levis", &fOutputData.levis);
  AddBranch("llen", &fOutputData.llen);
  AddBranch("lmomentum", &fOutputData.lmomentum);
  AddBranch("eccqe", &fOutputData.eccqe);

  AddBranch("track_pdg", &fOutputData.track_pdg);
  AddBranch("track_pdg_true", &fOutputData.track_pdg_true);
  AddBranch("track_evis", &fOutputData.track_evis);
  AddBranch("track_momentum", &fOutputData.track_momentum);
}


void LEESelection::Finalize() { 
  // Print out statistics
  // 
/*
  std::cout << "1e true: " << fCounters[kAny].true_1e
            << ", good: " << fCounters[kAny].good_1e
            << ", miss: " << fCounters[kAny].miss_1e
            << std::endl;

  std::cout << "1e eff: "
            << 1.0 * (fCounters[kAny].good_1e + fCounters[kAny].miss_1e) / fCounters[kAny].true_1e
            << std::endl;

  std::cout << "1e pur: "
            << 1.0 * fCounters[kAny].good_1e / (fCounters[kAny].good_1e + fCounters[kAny].miss_1e)
            << std::endl;

  std::cout << "1m true: " << fCounters[kAny].true_1m
            << ", good: " << fCounters[kAny].good_1m
            << ", miss: " << fCounters[kAny].miss_1m << std::endl;

  std::cout << "1m eff: "
            << 1.0 * (fCounters[kAny].good_1m + fCounters[kAny].miss_1m) / fCounters[kAny].true_1m
            << std::endl;

  std::cout << "1m pur: "
            << 1.0 * fCounters[kAny].good_1m / (fCounters[kAny].good_1m + fCounters[kAny].miss_1m)
            << std::endl;
*/
/*
  std::cout << "SHOWER, TRACK ENERGY RESOLUTION: "
            << _shower_energy_resolution << " "<< _track_energy_resolution << std::endl;

  std::cout << "ACCEPT NP " << _accept_np << std::endl;
  std::cout << "ACCEPT NTRK " << _accept_ntrk << std::endl;
*/
}

float LEESelection::nextTrackEnergyDistortion(float this_energy) {
  return this_energy * gRandom->Gaus(0, fConfig.track_energy_distortion);
}

float LEESelection::nextShowerEnergyDistortion(float this_energy) {
  return this_energy * gRandom->Gaus(0, fConfig.shower_energy_distortion);
}

int LEESelection::nextParticleID(float energy, int true_pdgid) {
  if (!fConfig.pdgid_matrix.is_set()) {
    return true_pdgid;
  }
  else {
    return fConfig.pdgid_matrix.get(energy)->particle_id(true_pdgid, gRandom->Uniform());
  }
}

bool LEESelection::ProcessEvent(gallery::Event& ev) {
  // Get handles for event data
  art::InputTag gtruth_tag(fConfig.mctruth_producer);

  auto const& gtruth_list = \
    *(ev.getValidHandle<std::vector<simb::GTruth> >(gtruth_tag));

  std::vector<evwgh::MCEventWeight> eventweights_list;
  if (fConfig.event_weight_producer != "") { 
    art::InputTag eventweight_tag(fConfig.event_weight_producer);
    eventweights_list = \
      *(ev.getValidHandle<std::vector<evwgh::MCEventWeight> >(eventweight_tag));
  }

  art::InputTag mctruth_tag(fConfig.mctruth_producer);
  auto const& mctruth_list = \
    *(ev.getValidHandle<std::vector<simb::MCTruth> >(mctruth_tag));

  art::InputTag mcshower_tag(fConfig.mcshower_producer);
  auto const& mcshower_list = \
    *(ev.getValidHandle<std::vector<sim::MCShower> >(mcshower_tag));

  art::InputTag mctrack_tag(fConfig.mctrack_producer);
  auto const& mctrack_list = \
    *(ev.getValidHandle<std::vector<sim::MCTrack> >(mctrack_tag));

  // Consistency checks
  assert(mctruth_list.size() == gtruth_list.size());

  // BNB flux weight
  double wbnb = 1.0;
  if (!eventweights_list.empty() &&
      eventweights_list[0].fWeight.find("bnbcorrection_FluxHist") != eventweights_list[0].fWeight.end()) {
    wbnb = eventweights_list[0].fWeight.at("bnbcorrection_FluxHist")[0];
  }

  // Loop through MC truth interactions
  if (mctruth_list.empty()) {
    return false;
  }

  // FIXME: Write a vector of OutputData objects to the output tree
  if (mctruth_list.size() > 1) {
    std::cerr << "LEESelection: Multiple neutrinos in event, processing first one only!" << std::endl;
  }

  const simb::MCTruth& mctruth = mctruth_list.at(0);
  const simb::GTruth& gtruth = gtruth_list.at(0);

  size_t ntracks = 0;
  size_t nshowers = 0;

  // Keep track of event particle content (currently a little redundant)
  std::vector<PIDParticle> particles_found;
  std::vector<PIDParticle> particles_true;

  // Loop through MC tracks
  for (size_t j=0; j<mctrack_list.size(); j++) {
    const sim::MCTrack& mct = mctrack_list.at(j);

    // Track properties
    double tlen = (mct.End().Position().Vect() - mct.Start().Position().Vect()).Mag();
    bool isEmpty = mct.empty();
    bool isFromNuVtx = util::IsFromNuVertex(mctruth, mct);
    bool isPrimary = mct.Process() == "primary";
    int tpdg = mct.PdgCode();
    float tke = mct.Start().E() - util::GetPDGMass(mct.PdgCode());

    // Track PID
    float energy_distortion = nextTrackEnergyDistortion(tke);
    int tpid = nextParticleID(tke + energy_distortion, tpdg);

    // Truth track cuts
    if (GoodObject(isFromNuVtx, isPrimary, tpdg, tke)) {
      particles_true.push_back({
        tpdg,
        tpdg,
        mct.Start().Momentum(),
        tke,
        util::ECCQE(mct.Start().Momentum()),
        tlen,
        !util::InFV(mct),
        mct.TrackID()
      });

      // Truth info on # of tracks
      ntracks++;
    }

    // PID track cuts
    if (GoodObject(isFromNuVtx, isPrimary, tpdg, tke + energy_distortion)) {
      particles_found.push_back({
        tpid,
        tpdg,
        mct.Start().Momentum(),
        tke + energy_distortion,
        util::ECCQE(mct.Start().Momentum(), energy_distortion),
        tlen,
        !util::InFV(mct),
        mct.TrackID()
      });
    }
  }

  // Loop through MC showers
  for (size_t j=0; j<mcshower_list.size(); j++) {
    const sim::MCShower& mcs = mcshower_list.at(j);

    // Shower properties
    double slen = (mcs.End().Position().Vect() - mcs.Start().Position().Vect()).Mag();
    bool isFromNuVtx = util::IsFromNuVertex(mctruth, mcs);
    bool isPrimary = mcs.Process() == "primary";
    int spdg = mcs.PdgCode();
    float ske = mcs.Start().E() - util::GetPDGMass(mcs.PdgCode());

    // Shower PID
    float energy_distortion = nextShowerEnergyDistortion(ske);
    int spid = nextParticleID(ske + energy_distortion, spdg);

    // Apply shower cuts
    if (GoodObject(isFromNuVtx, isPrimary, spdg, ske)) {
      particles_true.push_back({
        spdg,
        spdg,
        mcs.Start().Momentum(),
        ske,
        util::ECCQE(mcs.Start().Momentum()),
        slen,
        !util::InFV(mcs),
        mcs.TrackID()
      });

      // Truth info on # of showers
      nshowers++;
    }

    // PID shower cuts
    if (GoodObject(isFromNuVtx, isPrimary, spdg, ske + energy_distortion)) {
      particles_found.push_back({
        spid,
        spdg,
        mcs.Start().Momentum(),
        ske + energy_distortion,
        util::ECCQE(mcs.Start().Momentum(), energy_distortion),
        slen,
        !util::InFV(mcs),
        mcs.TrackID()
      });
    }
  }

  // Classify the event (found/true 1e/1m)
  bool f_1e = false;
  bool t_1e = false;
  bool f_1m = false;
  bool t_1m = false;

  for (auto sel : fConfig.selections) {
    f_1e |= TestSelection(particles_found, 11, sel);
    t_1e |= TestSelection(particles_true, 11, sel);
    f_1m |= TestSelection(particles_found, 13, sel);
    t_1m |= TestSelection(particles_true, 13, sel);
  }

  // Print out details for mis-IDs
  if ((f_1e && !t_1e) || (f_1m && !t_1m)) {
    std::cout << "true: " << mctruth.GetNeutrino().Nu().E() * 1000
              << "[" << mctruth.GetNeutrino().InteractionType() << "] ";
    for (size_t k=0; k<particles_true.size(); k++) {
      std::cout << particles_true[k] << " ";
    }
    std::cout << std::endl;
    std::cout << "pid: " << ntracks << " tracks, "
                         << nshowers << " showers; ";
    for (size_t k=0; k<particles_found.size(); k++) {
      std::cout << particles_found[k] << " ";
    }
    std::cout << std::endl;
  }

  // Fill output data
  fOutputData.np = GetNp(particles_found);
  fOutputData.ntrk = GetNtracks(particles_found);
  fOutputData.bnbweight = wbnb;
  fOutputData.dataset = fConfig.dataset_id;

  fOutputData.track_pdg.clear();
  fOutputData.track_pdg_true.clear();
  fOutputData.track_evis.clear();
  fOutputData.track_momentum.clear();

  for (size_t k=0; k<particles_found.size(); k++) {
    PIDParticle& p = particles_found[k];

    if (p.pdg == 11 || p.pdg == 13) {
      // Primary lepton
      fOutputData.lpid = p.pdg;
      fOutputData.lpdg = p.pdgtrue;
      fOutputData.lexit = p.exiting;
      fOutputData.levis = p.evis;
      fOutputData.lmomentum = p.p;
      fOutputData.eccqe = p.eccqe;
    }
    else {
      // Tracks
      fOutputData.track_pdg.push_back(p.pdg);
      fOutputData.track_pdg_true.push_back(p.pdgtrue);
      fOutputData.track_evis.push_back(p.evis);
      fOutputData.track_momentum.push_back(p.p);
    }
  }

  return f_1e || f_1m;
}

  }  // namespace LEETruthSelection
}  // namespace ana

// This line must be included for all selections!
DECLARE_SBN_PROCESSOR(ana::LEETruthSelection::LEESelection)

