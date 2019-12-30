#include "Profile.h"
#include "TH2D.h"
#include "TH3D.h"

namespace ana {
 namespace SBNOsc {

void TrackProfiles::Initialize(const std::string &postfix, unsigned nbinsx, double xlo, double xhi) {
#define TRACK_PROFILE(name, nbinsy, ylo, yhi) name = TH2Shared(new TH2D((#name "_" + postfix).c_str(), #name, nbinsx, xlo, xhi, nbinsy, ylo, yhi)); fAllHistos.push_back(name.Get());
#define TRACK_PROFILE3D(name, nbinsy, ylo, yhi, nbinsz, zlo, zhi) name = TH3Shared(new TH3D((#name "_" + postfix).c_str(), #name, nbinsx, xlo, xhi, nbinsy, ylo, yhi, nbinsz, zlo, zhi)); fAllHistos.push_back(name.Get())
  TRACK_PROFILE3D(range_v_true_mom, 50, 0., 2.5, 50, 0., 2.5);
  TRACK_PROFILE3D(mcs_v_true_mom, 50, 0., 2.5, 50, 0., 2.5);

  TRACK_PROFILE3D(range_minus_true, 50, 0., 2.5, 50, -1., 1.); 
  TRACK_PROFILE3D(mcs_minus_true, 50, 0., 2.5, 50, -1., 1.); 
  TRACK_PROFILE3D(pid_confusion_tr, 2, -0.5, 1.5, 2, -0.5, 1.5);
#undef TRACK_PROFILE3D
#undef TRACK_PROFILE
}


void TrackProfiles::Fill(const numu::ROOTValue &rootval, const numu::RecoTrack &track, const numu::RecoEvent &event) {
#define FILL3D(hist, x, y, z) hist.Fill(x, y, z)
  numu::LiteralValue literal = numu::MakeROOTLiteral(rootval, track, event);
  if (!literal.is_valid) return;

  float val = literal.is_float ? literal.data_num : (float) literal.data_int;

  if (track.match.has_match) {
    const numu::RecoTrack &true_track = event.true_tracks.at(track.match.mcparticle_id);

    FILL3D(range_v_true_mom, val, true_track.momentum, track.range_momentum);
    FILL3D(mcs_v_true_mom, val, true_track.momentum, track.mcs_momentum);

    FILL3D(range_minus_true, val, true_track.momentum, (track.range_momentum - true_track.momentum) / true_track.momentum);
    FILL3D(mcs_minus_true, val, true_track.momentum, (track.mcs_momentum - true_track.momentum) / true_track.momentum);
  }
  if (track.min_chi2 > 0) {
    bool is_proton_reco = track.chi2_proton < track.chi2_muon;
    if (track.match.has_match) {
      bool is_proton_true = abs(track.match.match_pdg) == 2212;
      bool is_muon_true = abs(track.match.match_pdg) == 13;
      if (is_proton_true || is_muon_true) {
        FILL3D(pid_confusion_tr, val, is_proton_true, is_proton_reco);
      }
    }
  }
#undef FILL3D
}

  } // end namespace SBNOsc
} // end namespace ana

