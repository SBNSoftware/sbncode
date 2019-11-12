#include "CosmicHisto.h"

#include "TH1D.h"
#include "TH2D.h"

namespace ana {
  namespace SBNOsc {

void CosmicHistos::Initialize(const std::string &postfix, const geo::BoxBoundedGeo &detector_volume) {
#define COSMIC_HISTO(name, n_bins, lo, hi) name = new TH1D((#name"_" + postfix).c_str(), #name, n_bins, lo, hi); fAllHistos.push_back(name)

  COSMIC_HISTO(enter_time, 100, -4000., 3000.);
  COSMIC_HISTO(enter_x, 400, detector_volume.MinX(), detector_volume.MaxX());
  COSMIC_HISTO(enter_y, 400, detector_volume.MinY(), detector_volume.MaxY());
  COSMIC_HISTO(enter_z, 500, detector_volume.MinZ(), detector_volume.MaxZ());
  COSMIC_HISTO(momentum, 100, 0., 5.);

#undef COSMIC_HISTO
}

void CosmicHistos::Fill(const std::vector<size_t> &cosmic_tracks, const std::map<size_t, numu::RecoTrack> &true_tracks) { 
  for (size_t id: cosmic_tracks) {
    const numu::RecoTrack &track = true_tracks.at(id);
    if (abs(track.pdgid) == 13) {
      enter_time->Fill(track.start_time);
      enter_x->Fill(track.start.X());
      enter_y->Fill(track.start.Y());
      enter_z->Fill(track.start.Z());
      momentum->Fill(track.momentum); 
    }
  }
}

  } // namespace SBNOsc
} // namespace ana

