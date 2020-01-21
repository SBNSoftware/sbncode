#include "CosmicHisto.h"

#include "TH1D.h"
#include "TH2D.h"

namespace ana {
  namespace SBNOsc {

void CosmicHistos::Initialize(const std::string &postfix, const geo::BoxBoundedGeo &detector_volume) {
#define COSMIC_HISTO(name, n_bins, lo, hi) name = new TH1D((#name + postfix).c_str(), #name, n_bins, lo, hi); StoreHisto(name)

  COSMIC_HISTO(enter_time, 600, -3000., 3000.);
  COSMIC_HISTO(enter_time_zoom, 300, -20., 10.);
  COSMIC_HISTO(enter_x, 400, detector_volume.MinX(), detector_volume.MaxX());
  COSMIC_HISTO(enter_y, 400, detector_volume.MinY(), detector_volume.MaxY());
  COSMIC_HISTO(enter_z, 500, detector_volume.MinZ(), detector_volume.MaxZ());
  COSMIC_HISTO(momentum, 100, 0., 5.);

#undef COSMIC_HISTO
}

void CosmicHistos::Fill(const std::map<size_t, numu::TrueParticle> &true_particles) { 
  for (const auto &pair: true_particles) {
    const numu::TrueParticle &particle = pair.second;
    // only plot cosmic muons
    if (particle.is_cosmic && abs(particle.pdgid) == 13) {
      enter_time->Fill(particle.start_time);
      enter_time_zoom->Fill(particle.start_time);
      enter_x->Fill(particle.start.X());
      enter_y->Fill(particle.start.Y());
      enter_z->Fill(particle.start.Z());
      momentum->Fill(particle.start_momentum.Mag()); 
    }
  }
}

  } // namespace SBNOsc
} // namespace ana

