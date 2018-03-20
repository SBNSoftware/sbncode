#include <algorithm>
#include <TBranch.h>
#include <TFile.h>
#include <TTree.h>
#include "gallery/Handle.h"
#include "canvas/Utilities/InputTag.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "json/json.h"
#include "Event.hh"
#include "Loader.hh"
#include "ProcessorBase.hh"

namespace core {

ProcessorBase::ProcessorBase() : fEventIndex(0), fOutputFilename("output.root") {}


ProcessorBase::~ProcessorBase() {}


void ProcessorBase::FillTree() {
  fTree->Fill();
  fEventIndex++;
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
    fOutputFilename = config->get("OutputFile", "output.root").asString();
  }

  // Open the output file and create the standard event tree
  fOutputFile = TFile::Open(fOutputFilename.c_str(), "recreate");
  fTree = new TTree("sbnana", "SBN Analysis Tree");
  fTree->AutoSave("overwrite");
  fEvent = new Event();
  fTree->Branch("events", &fEvent);
}


void ProcessorBase::Teardown() {
  // Write the standard tree and close the output file
  fOutputFile->cd();
  fTree->Write("sbnana", TObject::kOverwrite);
  fOutputFile->Close();
}


void ProcessorBase::BuildEventTree(gallery::Event& ev) {
  // Get MCTruth information
  auto const& mctruths = *ev.getValidHandle<std::vector<simb::MCTruth> >(fTruthTag);

  // Get MCEventWeight information
  gallery::Handle<std::vector<std::map<std::string,std::vector<double> > > > wgh;
  bool hasWeights = ev.getByLabel(fWeightTag, wgh);

  if (hasWeights) {
    assert(wgh->size() == mctruths.size());
  }

  fTree->GetEntry(fEventIndex);

  fEvent->ninteractions = std::min(Event::kMaxInteractions, mctruths.size());

  if (mctruths.size() > Event::kMaxInteractions) {
    std::cerr << "Too many interactions in event (" << mctruths.size() << "), "
              << "keeping maximum number (" << Event::kMaxInteractions << ")"
              << std::endl;
  }

  // Populate event tree
  for (size_t i=0; i<fEvent->ninteractions; i++) {
    Event::Interaction& interaction = fEvent->interactions[i];
    auto const& mctruth = mctruths.at(i);

    // Weights
    if (hasWeights) {
      interaction.weights = wgh->at(i);
    }

    // Neutrino
    const simb::MCNeutrino& nu = mctruth.GetNeutrino();
    interaction.neutrino.ccnc = nu.CCNC();
    interaction.neutrino.pdg = nu.Nu().PdgCode();
    interaction.neutrino.targetPDG = nu.Target();
    interaction.neutrino.intcode = nu.Mode();
    interaction.neutrino.bjorkenX = nu.X();
    interaction.neutrino.inelasticityY = nu.Y();
    interaction.neutrino.q2 = nu.QSqr();
    interaction.neutrino.w = nu.W();
    interaction.neutrino.energy = nu.Nu().EndMomentum().Energy();
    interaction.neutrino.momentum = nu.Nu().EndMomentum().Vect();

    // Primary lepton
    const simb::MCParticle& lepton = nu.Lepton();
    interaction.lepton.pdg = lepton.PdgCode();
    interaction.lepton.energy = lepton.Momentum(0).Energy();
    interaction.lepton.momentum = lepton.Momentum(0).Vect();

    // Final state system
    interaction.nfinalstate = \
      std::min((int)Event::kMaxFinalState, mctruth.NParticles());

    if (mctruth.NParticles() > Event::kMaxFinalState) {
      std::cerr << "Too many final state particles (" << mctruth.NParticles()
                << "), "
                << "keeping maximum number (" << Event::kMaxFinalState << ")"
                << std::endl;
    }

    for (int iparticle=0; iparticle<interaction.nfinalstate; iparticle++) {
      Event::FinalStateParticle& fsp = interaction.finalstate[iparticle];
      const simb::MCParticle& particle = mctruth.GetParticle(iparticle);

      if (particle.Process() != "primary") {
        continue;
      }

      fsp.pdg = particle.PdgCode();
      fsp.energy = particle.Momentum(0).Energy();
      fsp.momentum = particle.Momentum(0).Vect();
    }
  }
}

}  // namespace core

