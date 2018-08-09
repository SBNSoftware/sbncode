#include <iostream>
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "core/Event.hh"

namespace ana {
  namespace SBNOsc {

void hello() {
  std::cout << "Hello SBNOsc!" << std::endl;
}


Event::Interaction TruthReco(const simb::MCTruth& mctruth) {
  Event::Interaction interaction;

  // Neutrino
  const simb::MCNeutrino& nu = mctruth.GetNeutrino();
  interaction.neutrino.energy = nu.Nu().EndMomentum().Energy();

  // Primary lepton
  const simb::MCParticle& lepton = nu.Lepton();
  interaction.lepton.pdg = lepton.PdgCode();
  interaction.lepton.energy = lepton.Momentum(0).Energy();
  interaction.lepton.momentum = lepton.Momentum(0).Vect();

  // Hadronic system
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

  return interaction;
}

  }  // namespace SBNOsc
}  // namespace ana

