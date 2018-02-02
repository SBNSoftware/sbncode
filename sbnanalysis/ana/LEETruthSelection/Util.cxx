#include <string>
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/MCBase/MCShower.h"
#include "lardataobj/MCBase/MCTrack.h"
#include "TH2F.h"
#include "TDatabasePDG.h"
#include "Util.h"

namespace ana {
  namespace LEETruthSelection {
    namespace util {

// From J. Zennamo's pion selection
const double fv_x_lo =    0.0, fv_x_hi =  256.35;
const double fv_y_lo = -116.5, fv_y_hi =  116.50;
const double fv_z_lo =    0.0, fv_z_hi = 1036.80;


bool InFV(const sim::MCShower& s) {
  return (s.End().X() > fv_x_lo && s.End().X() < fv_x_hi &&
          s.End().Y() > fv_y_lo && s.End().Y() < fv_y_hi &&
          s.End().Z() > fv_z_lo && s.End().Z() < fv_z_hi);
}


bool InFV(const sim::MCTrack& t) {
  return (t.End().X() > fv_x_lo && t.End().X() < fv_x_hi &&
          t.End().Y() > fv_y_lo && t.End().Y() < fv_y_hi &&
          t.End().Z() > fv_z_lo && t.End().Z() < fv_z_hi);
}


bool IsFromNuVertex(const simb::MCTruth& mc, const sim::MCShower& show,
                    float distance) {
  TLorentzVector nuVtx = mc.GetNeutrino().Nu().Trajectory().Position(0);
  TLorentzVector showStart = show.Start().Position();
  return (showStart - nuVtx).Mag() < distance;
}


bool IsFromNuVertex(const simb::MCTruth& mc, const sim::MCTrack& track,
                    float distance) {
  TLorentzVector nuVtx = mc.GetNeutrino().Nu().Trajectory().Position(0);
  TLorentzVector trkStart = track.Start().Position();
  return (trkStart - nuVtx).Mag() < distance;
}


double ECCQE(const TLorentzVector& inp_v, float energy_distortion) {
  double e_dist = (double) energy_distortion;
  // Based on D. Kaleko, LowEnergyExcess LArLite module ECCQECalculator
  double M_n = 939.565; // MeV/c^2
  double M_p = 938.272; // MeV/c^2
  double M_e = 0.511; // MeV/c^2
  double bindingE = 30.0; // MeV

  double mp2 = M_p * M_p;
  double me2 = M_e * M_e;
  double mnb = M_n - bindingE;

  TLorentzVector v(inp_v);

  double l_energy = v.E() + e_dist;

  double l_mom = sqrt(l_energy*l_energy - me2);

  double l_theta = \
    acos(v.Pz() / sqrt(v.Px()*v.Px() + v.Py()*v.Py() + v.Pz()*v.Pz()));

  double enu_top = mp2 - mnb*mnb - me2 + 2.0 * mnb * l_energy;
  double enu_bot = 2.0 * (mnb - l_energy + l_mom * cos(l_theta));

  return enu_top / enu_bot;
}


double GetPDGMass(const int pdg) {
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

    }  // namespace util
  }  // namespace LEETruthSelection
}  // namespace ana

