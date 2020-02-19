#include "TrueParticleHisto.h"

#include "TH1D.h"
#include "TH2D.h"

#include "TDatabasePDG.h"

namespace ana {
  namespace SBNOsc {

void TrueParticleHistos::Initialize(const std::string &postfix) {
#define PARTICLE_HISTO(name, n_bins, lo, hi) name = new TH1D((#name + postfix).c_str(), #name, n_bins, lo, hi); StoreHisto(name)

  PARTICLE_HISTO(momentum, 100, 0, 4);
  PARTICLE_HISTO(deposited_energy, 100, 0, 2);
  PARTICLE_HISTO(length, 100, 0, 600);
  PARTICLE_HISTO(kinetic_energy, 100, 0, 2);
  PARTICLE_HISTO(frac_deposited, 100, 0., 3.);

#undef PARTICLE_HISTO
}

double PDGMass(int pdg) {
  static const TDatabasePDG *PDGTable(new TDatabasePDG);

  // regular particle
  if (pdg < 1000000000) {
    TParticlePDG* ple = PDGTable->GetParticle(pdg);
    if (ple == NULL) return -1;
    return ple->Mass();
  }
  // ion
  else {
    int p = (pdg % 10000000) / 10000;
    int n = (pdg % 10000) / 10 - p;
    return (PDGTable->GetParticle(2212)->Mass() * p +
            PDGTable->GetParticle(2112)->Mass() * n);
  }
}

void TrueParticleHistos::Fill(const numu::TrueParticle &particle) { 
  momentum->Fill(particle.start_momentum.Mag());
  deposited_energy->Fill(particle.deposited_energy);
  length->Fill(particle.length);

  float mass = PDGMass(particle.pdgid);
  if (mass > 0) {
    float p = particle.start_momentum.Mag();
    float ke = sqrt(mass*mass + p*p) - mass;
    kinetic_energy->Fill(ke);
    frac_deposited->Fill(particle.deposited_energy / ke);
  }
}

  } // namespace SBNOsc
} // namespace ana

