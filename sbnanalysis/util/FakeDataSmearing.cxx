#include "core/Event.hh"
#include "core/Experiment.hh"
#include "FakeDataSmearing.hh"

namespace util {

size_t CCNPi(const event::Interaction* truth) {
  if (!truth->neutrino.iscc) {
    return false;
  }

  size_t npi = 0;
  for (size_t i=0; i<truth->finalstate.size(); i++) {
    if (truth->finalstate[i].status_code == 1 &&
        (abs(truth->finalstate[i].pdg) == 211 || truth->finalstate[i].pdg == 11)) {
      npi++;
    }
  }

  return npi;
}


void ApplyFakeDataSmearing(SmearingMode mode,
                           Experiment exp,
                           const event::Interaction* truth,
                           event::RecoInteraction* reco) {
  for (int i=kMode1; i<=mode; i++) {
    // Dataset 1
    if (i == kMode1) {
      // Do nothing
    }

    // Dataset 2
    else if (i == kMode2) {
      // Decrease CCQE by 10% for Q^2 < 0.3 GeV^2 and decrease by 3% in remaining Q^2 range
      if (truth->neutrino.iscc && truth->neutrino.genie_intcode == 0) {
        if (truth->neutrino.Q2 < 0.3) {
          reco->weight *= 0.9;
        }
        else {
          reco->weight *= 0.97;
        }
      }

      // Increase CCMEC by 50%
      if (truth->neutrino.iscc && truth->neutrino.genie_intcode == 10) {
        reco->weight *= 1.5;
      }

      // Throw away a random 20% of CC1pi events where the pi is subsequently absorbed
      // For now, throw away 20% of pionless true CCRes
      if (truth->neutrino.iscc && truth->neutrino.genie_intcode == 1 && CCNPi(truth) == 0) {
        reco->weight *= 0.8;
      }

      // Increase NC bkg by 30%
      if (truth->neutrino.isnc) {
        reco->weight *= 1.3;
      }
    }

    // Dataset 3
    else if (i == kMode3) {
      // Decrease muon-neutrino flux below 1 GeV by 2% at SBND ...
      if (exp == kExpSBND && truth->neutrino.pdg == 14) {
        reco->weight *= 0.98;
      }

      // ... and by 1% at ICARUS
      if (exp == kExpICARUS && truth->neutrino.pdg == 14) {
        reco->weight *= 0.99;
      }
    }

    // Dataset 4
    else if (i == kMode4) {
      float ccqe_eshift = 0;
      float ccmec_eshift = 0;
      float cc1pi_eshift = 0;

      if (exp == kExpSBND) {
        float ccqe_eshift = -0.02;
        float ccmec_eshift = -0.03;
        float cc1pi_eshift = -0.04;
      }
      else if (exp == kExpICARUS) {
        float ccqe_eshift = -0.01;
        float ccmec_eshift = -0.02;
        float cc1pi_eshift = -0.03;
      }

      // In the calculation of Ereco, remove x MeV from CCQE events ...
      if (truth->neutrino.iscc && truth->neutrino.genie_intcode == 0) {
        reco->reco_energy += ccqe_eshift;
      }

      // CCMEC events
      if (truth->neutrino.iscc && truth->neutrino.genie_intcode == 10) {
        reco->reco_energy += ccmec_eshift;
      }

      // CC1pi events
      if (CCNPi(truth) == 1) {
        reco->reco_energy += cc1pi_eshift;
      }
    }

    // Dataset 5
    else if (i == kMode5) {
      // Apply muon neutrino disappearance
      if (truth->neutrino.pdg == 14) {
        float sinsq2tmm = 0.08;  // sin^2(2theta_{mumu})
        float dm41sq = 2;  // Delta m^2_{41} [eV^2]
        float l = truth->neutrino.baseline / 1000;  // [km]
        float enu = truth->neutrino.energy;  // [GeV]
        float pmm = 1.0 - sinsq2tmm * pow(sin(1.27 * dm41sq * l / enu), 2);

        reco->weight *= pmm;
      }
    }
  }
}

}  // namespace util

