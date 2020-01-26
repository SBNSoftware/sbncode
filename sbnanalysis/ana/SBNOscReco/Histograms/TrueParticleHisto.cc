#include "TrueParticleHisto.h"

#include "TH1D.h"
#include "TH2D.h"

namespace ana {
  namespace SBNOsc {

void TrueParticleHistos::Initialize(const std::string &postfix) {
#define PARTICLE_HISTO(name, n_bins, lo, hi) name = new TH1D((#name + postfix).c_str(), #name, n_bins, lo, hi); StoreHisto(name)

  PARTICLE_HISTO(momentum, 100, 0, 4);
  PARTICLE_HISTO(deposited_energy, 100, 0, 2);
  PARTICLE_HISTO(length, 100, 0, 600);

#undef PARTICLE_HISTO
}

void TrueParticleHistos::Fill(const numu::TrueParticle &particle) { 
  momentum->Fill(particle.start_momentum.Mag());
  deposited_energy->Fill(particle.deposited_energy);
  length->Fill(particle.length);
}

  } // namespace SBNOsc
} // namespace ana

