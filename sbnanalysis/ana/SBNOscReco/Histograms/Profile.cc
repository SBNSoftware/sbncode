#include "Profile.h"
#include "TH2D.h"

namespace ana {
 namespace SBNOsc {

void TrackProfiles::Initialize(const std::string &postfix, unsigned nbinsx, double xlo, double xhi) {
#define TRACK_PROFILE(name, nbinsy, ylo, yhi) name = new TH2D((#name "_" + postfix).c_str(), #name, nbinsx, xlo, xhi, nbinsy, ylo, yhi); all_histos.push_back(name)
  TRACK_PROFILE(range_bias, 50, -1., 1.);
  TRACK_PROFILE(range_var, 50., 0., 1.);
  TRACK_PROFILE(mcs_bias, 50, -1., 1.);
  TRACK_PROFILE(mcs_var, 50., 0., 1.);
#undef TRACK_PROFILE
}


void TrackProfiles::Fill(const numu::ROOTValue &rootval, const numu::RecoTrack &track, const numu::RecoEvent &event) {
  numu::LiteralValue literal = numu::MakeROOTLiteral(rootval, track, event);
  if (!literal.is_valid) return;

  float val = literal.is_float ? literal.data_num : (float) literal.data_int;

  if (track.match.has_match) {
    const numu::RecoTrack &true_track = event.true_tracks.at(track.match.mcparticle_id);

    range_bias->Fill(val, (track.range_momentum - true_track.momentum) / true_track.momentum);
    range_var->Fill(val, (track.range_momentum - true_track.momentum) * (track.range_momentum - true_track.momentum) / (true_track.momentum * true_track.momentum));

    mcs_bias->Fill(val, (track.fwd_mcs_momentum - true_track.momentum) / true_track.momentum);
    mcs_var->Fill(val, (track.fwd_mcs_momentum - true_track.momentum) * (track.fwd_mcs_momentum - true_track.momentum) / (true_track.momentum * true_track.momentum));
  }
}

  } // end namespace SBNOsc
} // end namespace ana

