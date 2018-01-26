#include <TFile.h>
#include <TTree.h>
#include "gallery/ValidHandle.h"
#include "canvas/Utilities/InputTag.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "Event.hh"
#include "ProcessorBase.hh"

namespace io {

ProcessorBase::ProcessorBase() : fEventIndex(0), fOutputFilename("output.root") {}


ProcessorBase::~ProcessorBase() {}


void ProcessorBase::Initialize() {
  fOutputFile = TFile::Open(fOutputFilename.c_str(), "recreate");
  fTree = new TTree("sbnana", "SBN Analysis Tree");
  fEvent = new Event();
  fTree->Branch("events", &fEvent);
}


void ProcessorBase::Finalize() {
  fOutputFile->cd();
  fTree->Write();
  fOutputFile->Close();
}


void ProcessorBase::ProcessFile(std::vector<std::string> filenames) {
  for (gallery::Event ev(filenames); !ev.atEnd(); ev.next()) {
    FillEventTree(ev);
    ProcessEvent(ev);
    fTree->Fill();
    fEventIndex++;
  }
}


void ProcessorBase::FillEventTree(gallery::Event& ev) {
  art::InputTag mctruths_tag = { "generator" };
 
  auto const& mctruths = *ev.getValidHandle<std::vector<simb::MCTruth> >(mctruths_tag);

  fTree->GetEntry(fEventIndex);

  // Populate event tree
  for (size_t i=0; i<mctruths.size(); i++) {
    Event::Interaction interaction;
    auto const& mctruth = mctruths.at(i);

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

    // Hadronic system
    for (int iparticle=0; iparticle<mctruth.NParticles(); iparticle++) {
      const simb::MCParticle& particle = mctruth.GetParticle(iparticle);

      if (particle.Process() != "primary") {
        continue;
      }

      Event::FinalStateParticle fsp;
      fsp.pdg = particle.PdgCode();
      fsp.energy = particle.Momentum(0).Energy();
      fsp.momentum = particle.Momentum(0).Vect();

      interaction.hadrons.push_back(fsp);
    }

    fEvent->interactions.push_back(interaction);
  }
}

}  // namespace io

