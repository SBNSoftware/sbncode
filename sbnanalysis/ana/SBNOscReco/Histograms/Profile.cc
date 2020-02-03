#include "Profile.h"
#include "TH2D.h"
#include "TH3D.h"

#include "../NumuReco/TrackAlgo.h"

namespace ana {
 namespace SBNOsc {

void TrackProfiles::Initialize(const std::string &postfix, unsigned nbinsx, double xlo, double xhi) {
#define TRACK_PROFILE(name, nbinsy, ylo, yhi) name = new TH2D((#name  + postfix).c_str(), #name, nbinsx, xlo, xhi, nbinsy, ylo, yhi); StoreHisto(name);
#define TRACK_PROFILE3D(name, nbinsy, ylo, yhi, nbinsz, zlo, zhi) name = new TH3D((#name  + postfix).c_str(), #name, nbinsx, xlo, xhi, nbinsy, ylo, yhi, nbinsz, zlo, zhi); StoreHisto(name)
  TRACK_PROFILE3D(range_v_true_mom, 50, 0., 2.5, 50, 0., 2.5);
  TRACK_PROFILE3D(mcs_v_true_mom, 50, 0., 2.5, 50, 0., 2.5);

  TRACK_PROFILE3D(range_minus_true, 50, 0., 2.5, 50, -1., 1.); 
  TRACK_PROFILE3D(mcs_minus_true, 50, 0., 2.5, 50, -1., 1.); 
  TRACK_PROFILE3D(pid_confusion_tr, 2, -0.5, 1.5, 2, -0.5, 1.5);
#undef TRACK_PROFILE3D
#undef TRACK_PROFILE
}


void TrackProfiles::Fill(float val, const numu::RecoTrack &track, const numu::RecoEvent &event) {
  int mcparticle_id = track.truth.GetPrimaryMatchID();

  if (mcparticle_id >= 0) {
    const numu::TrueParticle &true_particle = event.particles.at(mcparticle_id);

    range_v_true_mom->Fill(val, true_particle.start_momentum.Mag(), numu::RangeMomentum(track));
    mcs_v_true_mom->Fill(val, true_particle.start_momentum.Mag(), numu::MCSMomentum(track));

    range_minus_true->Fill(val, true_particle.start_momentum.Mag(), (numu::RangeMomentum(track) - true_particle.start_momentum.Mag()) / true_particle.start_momentum.Mag());
    mcs_minus_true->Fill(val, true_particle.start_momentum.Mag(), (numu::MCSMomentum(track) - true_particle.start_momentum.Mag()) / true_particle.start_momentum.Mag());
    if (track.min_chi2 > 0) {
      bool is_proton_reco = track.chi2_proton < track.chi2_muon;
      bool is_proton_true = abs(true_particle.pdgid) == 2212;
      bool is_muon_true = abs(true_particle.pdgid) == 13;
      if (is_proton_true || is_muon_true) {
        pid_confusion_tr->Fill(val, is_proton_true, is_proton_reco);
      }
    }
  }
}

  } // end namespace SBNOsc
} // end namespace ana

