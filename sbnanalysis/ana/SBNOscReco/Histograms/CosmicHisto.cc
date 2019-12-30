#include "CosmicHisto.h"

#include "TH1D.h"
#include "TH2D.h"

namespace ana {
  namespace SBNOsc {

void CosmicHistos::Initialize(const std::string &postfix, const geo::BoxBoundedGeo &detector_volume) {
#define COSMIC_HISTO(name, n_bins, lo, hi) name = TH1Shared(new TH1D((#name"_" + postfix).c_str(), #name, n_bins, lo, hi)); fAllHistos.push_back(name.Get())

  COSMIC_HISTO(enter_time, 600, -3000., 3000.);
  COSMIC_HISTO(enter_time_zoom, 300, -20., 10.);
  COSMIC_HISTO(enter_x, 400, detector_volume.MinX(), detector_volume.MaxX());
  COSMIC_HISTO(enter_y, 400, detector_volume.MinY(), detector_volume.MaxY());
  COSMIC_HISTO(enter_z, 500, detector_volume.MinZ(), detector_volume.MaxZ());
  COSMIC_HISTO(momentum, 100, 0., 5.);

#undef COSMIC_HISTO
}

void CosmicHistos::Fill(const std::vector<size_t> &cosmic_tracks, const std::map<size_t, numu::RecoTrack> &true_tracks) { 
#define FILL(hist, val) hist.Fill(val);
  for (size_t id: cosmic_tracks) {
    const numu::RecoTrack &track = true_tracks.at(id);
    if (abs(track.pdgid) == 13) {
      FILL(enter_time, track.start_time);
      FILL(enter_time_zoom, track.start_time);
      FILL(enter_x, track.start.X());
      FILL(enter_y, track.start.Y());
      FILL(enter_z, track.start.Z());
      FILL(momentum, track.momentum);
    }
  }
#undef FILL
}

  } // namespace SBNOsc
} // namespace ana

