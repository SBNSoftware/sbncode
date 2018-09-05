#include <algorithm>
#include <TBranch.h>
#include <TFile.h>
#include <TTree.h>
#include "gallery/ValidHandle.h"
#include "gallery/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "nusimdata/SimulationBase/GTruth.h"
#include "larsim/EventWeight/Base/MCEventWeight.h"
#include "json/json.h"
#include "Event.hh"
#include "Loader.hh"
#include "util/Interaction.hh"
#include "ProcessorBase.hh"

namespace core {

ProcessorBase::ProcessorBase() : fEventIndex(0), fOutputFilename("output.root") {}


ProcessorBase::~ProcessorBase() {}


void ProcessorBase::FillTree() {
  fTree->Fill();
  fEventIndex++;
}

void ProcessorBase::EventCleanup() {
  fEvent->metadata.Init();
  fEvent->truth.clear();
  fEvent->reco.clear();
}


void ProcessorBase::Initialize(char* config) {
  Json::Value* cfg = LoadConfig(config);
  Initialize(cfg);
}


void ProcessorBase::Setup(char* config) {
  Json::Value* cfg = LoadConfig(config);
  Setup(cfg);
}


void ProcessorBase::Setup(Json::Value* config) {
  // Load configuration parameters
  fTruthTag = { "generator" };
  fWeightTag = { "eventweight" };

  if (config) {
    fTruthTag = { config->get("MCTruthTag", "generator").asString() };
    fWeightTag = { config->get("MCWeightTag", "eventweight").asString() };
    fMCTrackTag = { config->get("MCTrackTag", "mcreco").asString() };
    fMCShowerTag = { config->get("MCShowerTag", "mcreco").asString() };
    fMCParticleTag = { config->get("MCParticleTag", "largeant").asString() };
    fOutputFilename = config->get("OutputFile", "output.root").asString();
  }

  // Open the output file and create the standard event tree
  fOutputFile = TFile::Open(fOutputFilename.c_str(), "recreate");
  fTree = new TTree("sbnana", "SBN Analysis Tree");
  fTree->AutoSave("overwrite");
  fEvent = new Event();
  fTree->Branch("events", &fEvent);
  fReco = &fEvent->reco;
}


void ProcessorBase::Teardown() {
  // Write the standard tree and close the output file
  fOutputFile->cd();
  fTree->Write("sbnana", TObject::kOverwrite);
  fOutputFile->Close();
}


void ProcessorBase::BuildEventTree(gallery::Event& ev) {
  // Get MCTruth information
  auto const& mctruths = \
    *ev.getValidHandle<std::vector<simb::MCTruth> >(fTruthTag);
  
  gallery::Handle<std::vector<simb::GTruth> > gtruths_handle;
  ev.getByLabel(fTruthTag,gtruths_handle);
  bool genie_truth_is_valid = gtruths_handle.isValid();

  // Get MCEventWeight information
  gallery::Handle<std::vector<evwgh::MCEventWeight> > wgh;
  bool hasWeights = ev.getByLabel(fWeightTag, wgh);

  if (hasWeights) {
    assert(wgh->size() == mctruths.size());
  }

  fTree->GetEntry(fEventIndex);

  // Populate event tree
  for (size_t i=0; i<mctruths.size(); i++) {
    Event::Interaction interaction;

    auto const& mctruth = mctruths.at(i);

    // TODO: What to do with cosmic MC?
    // For now, ignore them
    if (!mctruth.NeutrinoSet()) continue;

    // Weights
    if (hasWeights) {
      interaction.weights = wgh->at(i).fWeight;
    }

    // Neutrino
    const simb::MCNeutrino& nu = mctruth.GetNeutrino();
    interaction.neutrino.isnc =   nu.CCNC()  && (nu.Mode() != simb::kWeakMix);
    interaction.neutrino.iscc = (!nu.CCNC()) && (nu.Mode() != simb::kWeakMix);
    interaction.neutrino.pdg = nu.Nu().PdgCode();
    interaction.neutrino.targetPDG = nu.Target();
    interaction.neutrino.genie_intcode = nu.Mode();
    interaction.neutrino.bjorkenX = nu.X();
    interaction.neutrino.inelasticityY = nu.Y();
    interaction.neutrino.Q2 = nu.QSqr();
    interaction.neutrino.w = nu.W();
    interaction.neutrino.energy = nu.Nu().EndMomentum().Energy();
    interaction.neutrino.momentum = nu.Nu().EndMomentum().Vect();

    // Primary lepton
    const simb::MCParticle& lepton = nu.Lepton();
    interaction.lepton.pdg = lepton.PdgCode();
    interaction.lepton.energy = lepton.Momentum(0).Energy();
    interaction.lepton.momentum = lepton.Momentum(0).Vect();

    // get CCQE energy from lepton info
    interaction.neutrino.eccqe = util::ECCQE(interaction.lepton.momentum, interaction.lepton.energy);

    // Hadronic system
    TLorentzVector q_labframe = nu.Nu().EndMomentum() - lepton.Momentum(0);
    interaction.neutrino.q0_lab = q_labframe.E();
    interaction.neutrino.modq_lab = q_labframe.P();

    for (int iparticle=0; iparticle<interaction.finalstate.size(); iparticle++) {
      Event::FinalStateParticle fsp;
      const simb::MCParticle& particle = mctruth.GetParticle(iparticle);

      if (particle.Process() != "primary") {
        continue;
      }

      fsp.pdg = particle.PdgCode();
      fsp.energy = particle.Momentum(0).Energy();
      fsp.momentum = particle.Momentum(0).Vect();

      interaction.finalstate.push_back(fsp);
    }

    // GENIE specific
    if (genie_truth_is_valid) {
      auto const& gtruth = gtruths_handle->at(i);
      TLorentzVector q_nucframe(q_labframe);
      // This nucleon momentum should be added to MCNeutrino so we don't
      // have to rely on GTruth
      const TLorentzVector& nucP4 = gtruth.fHitNucP4;
      TVector3 nuc_boost(nucP4.BoostVector());
      q_nucframe.Boost(nuc_boost);
      interaction.neutrino.modq = q_nucframe.P();
      interaction.neutrino.q0 = q_nucframe.E();
    }

    fEvent->truth.push_back(interaction);
  }
}

}  // namespace core

