#include <string>
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "TH2F.h"
#include "TDatabasePDG.h"
#include "TSUtil.h"

namespace ana {
namespace lee_truth_selection {
namespace tsutil {

// From J. Zennamo's pion selection
const double fv_x_lo =    0.0, fv_x_hi =  256.35;
const double fv_y_lo = -116.5, fv_y_hi =  116.50;
const double fv_z_lo =    0.0, fv_z_hi = 1036.80;


bool inFV(const sim::MCShower& s) {
  return (s.End().X() > fv_x_lo && s.End().X() < fv_x_hi &&
          s.End().Y() > fv_y_lo && s.End().Y() < fv_y_hi &&
          s.End().Z() > fv_z_lo && s.End().Z() < fv_z_hi);
}


bool inFV(const sim::MCTrack& t) {
  return (t.End().X() > fv_x_lo && t.End().X() < fv_x_hi &&
          t.End().Y() > fv_y_lo && t.End().Y() < fv_y_hi &&
          t.End().Z() > fv_z_lo && t.End().Z() < fv_z_hi);
}


bool isFromNuVertex(const simb::MCTruth& mc, const sim::MCShower& show,
                    float distance) {
  TLorentzVector nuVtx = mc.GetNeutrino().Nu().Trajectory().Position(0);
  TLorentzVector showStart = show.Start().Position();
  return (showStart - nuVtx).Mag() < distance;
}


bool isFromNuVertex(const simb::MCTruth& mc, const sim::MCTrack& track,
                    float distance) {
  TLorentzVector nuVtx = mc.GetNeutrino().Nu().Trajectory().Position(0);
  TLorentzVector trkStart = track.Start().Position();
  return (trkStart - nuVtx).Mag() < distance;
}


double eccqe(const TLorentzVector& inp_v, float energy_distortion, float angle_distortion) {
  double e_dist = (double) energy_distortion;
  // Based on D. Kaleko, LowEnergyExcess LArLite module ECCQECalculator
  double M_n = 939.565; // MeV/c^2
  double M_p = 938.272; // MeV/c^2
  double M_e = 0.511; // MeV/c^2
  double bindingE = 30.0; // MeV

  double mp2 = M_p * M_p;
  double me2 = M_e * M_e;
  double mnb = M_n - bindingE;

  // mess around with lorentz vector
  TLorentzVector v(inp_v);
  v.SetTheta( v.Theta() + angle_distortion);

  double l_energy = v.E() + e_dist;

  double l_mom = sqrt(l_energy*l_energy - me2);

  double l_theta = \
    acos(v.Pz() / sqrt(v.Px()*v.Px() + v.Py()*v.Py() + v.Pz()*v.Pz()));

  double enu_top = mp2 - mnb*mnb - me2 + 2.0 * mnb * l_energy;
  double enu_bot = 2.0 * (mnb - l_energy + l_mom * cos(l_theta));

  return enu_top / enu_bot;
}


double get_pdg_mass(const int pdg) {
  if (!gPDGTable) {
    gPDGTable = new TDatabasePDG;
  }
  if (pdg < 1000000000) {
    // A regular old particle
    TParticlePDG* ple = gPDGTable->GetParticle(pdg);
    return ple->Mass() * 1000.0;
  }
  else {
    // Ion, total up proton and neutron masses
    int p = (pdg % 10000000) / 10000;
    int n = (pdg % 10000) / 10 - p;
    return (gPDGTable->GetParticle(2212)->Mass() * p +
            gPDGTable->GetParticle(2112)->Mass() * n) * 1000.0;
  }
}

bool is_shower_pdgid(int pdg) {
  return (pdg == 11 || pdg == 22);
}

TH2F* HistDEdx(const sim::MCTrack& t, const std::string name, int lowbin) {
  TH2F* h = new TH2F(name.c_str(), ";Residual range (cm?);dE/dx (MeV/cm)",
                     100, 0, 200, 100, 0, 10);

  double s = 0;  // Track length
  TLorentzVector pos = t.End().Position();  // Start from end, work back

  for (long k=t.size()-2; k>=0; k--) {
    double dedx = t.dEdx()[k];
    if (h->GetXaxis()->FindBin(s) >= lowbin &&
        h->GetYaxis()->FindBin(dedx) >= lowbin) {
      h->Fill(s, dedx);
    }

    s += (pos.Vect() - t[k].Position().Vect()).Mag();
    pos = t[k].Position();
  }

  return h;
}

}  // namespace tsutil
}
}
