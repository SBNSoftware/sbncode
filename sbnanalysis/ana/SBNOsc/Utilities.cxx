#include <iostream>

#include "TDatabasePDG.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "core/Event.hh"

#include "Utilities.h"

#include "ubcore/LLBasicTool/GeoAlgo/GeoAABox.h"
#include "ubcore/LLBasicTool/GeoAlgo/GeoAlgo.h"
#include "ubcore/LLBasicTool/GeoAlgo/GeoLineSegment.h"

namespace ana {
  namespace SBNOsc {

void hello() {
  std::cout << "Hello SBNOsc!" << std::endl;
}


Event::Interaction TruthReco(const simb::MCTruth& mctruth) {
  Event::Interaction interaction;

  // Neutrino
  const simb::MCNeutrino& nu = mctruth.GetNeutrino();
  interaction.neutrino.energy = nu.Nu().EndMomentum().Energy();

  // Primary lepton
  const simb::MCParticle& lepton = nu.Lepton();
  interaction.lepton.pdg = lepton.PdgCode();
  interaction.lepton.energy = lepton.Momentum(0).Energy();
  interaction.lepton.momentum = lepton.Momentum(0).Vect();

  // Hadronic system
  for (int iparticle=0; iparticle<interaction.finalstate.size(); iparticle++) {
    Event::FinalStateParticle fsp;
    const simb::MCParticle& particle = mctruth.GetParticle(iparticle);

    if (particle.Process() != "primary") {
      continue;
    }

    fsp.pdg = particle.PdgCode();
    fsp.energy = particle.Momentum(0).Energy();
    fsp.momentum = particle.Momentum(0).Vect();

    interaction.finalstate.push_back(fsp);
  }

  return interaction;
}

double ECCQE(const TVector3& l_momentum, double l_energy, double energy_distortion, double angle_distortion) {
  // Based on D. Kaleko, LowEnergyExcess LArLite module ECCQECalculator
  double M_n = 0.939565; // GeV/c^2
  double M_p = 0.938272; // GeV/c^2
  double M_e = 0.000511; // GeV/c^2
  double bindingE = 0.0300; // GeV
  
  double mp2 = M_p * M_p;
  double me2 = M_e * M_e;
  double mnb = M_n - bindingE;

  // mess around with lorentz vector
  TVector3 v(l_momentum);
  v.SetTheta( v.Theta() + angle_distortion);
  l_energy = l_energy + energy_distortion;
  double l_mom = sqrt(l_energy*l_energy - me2);
  double l_theta = \
    acos(v.Pz() / sqrt(v.Px()*v.Px() + v.Py()*v.Py() + v.Pz()*v.Pz()));
  double enu_top = mp2 - mnb*mnb - me2 + 2.0 * mnb * l_energy;
  double enu_bot = 2.0 * (mnb - l_energy + l_mom * cos(l_theta));
  return enu_top / enu_bot;
}

double NuMuOscillation(double numu_energy, double numu_dist, double osc_dm2, double osc_angle) {
  double overlap = sin(2*osc_angle);
  double energy_factor = sin(1.27 * osc_dm2 * numu_dist / numu_energy);
  return 1 - overlap * overlap * energy_factor * energy_factor;

}

double containedLength(const TVector3 &v0, const TVector3 &v1, const std::vector<geoalgo::AABox> &boxes) {
  static const geoalgo::GeoAlgo algo;

  // if points are the same, return 0
  if ((v0 - v1).Mag() < 1e-6) return 0;

  // construct individual points
  geoalgo::Point_t p0(v0);
  geoalgo::Point_t p1(v1);

  // construct line segment
  geoalgo::LineSegment line(p0, p1);

  double length = 0;

  // total contained length is sum of lengths in all boxes
  // assuming they are non-overlapping
  for (auto const &box: boxes) {
    int n_contained = box.Contain(p0) + box.Contain(p1);
    // both points contained -- length is total length (also can break out of loop)
    if (n_contained == 2) {
      length = (v1 - v0).Mag();
      break;
    }
    // one contained -- have to find intersection point (which must exist) 
    if (n_contained == 1) {
      auto intersections = algo.Intersection(line, box);
      assert(intersections.size() == 1); // must have one intersection
      // get TVector at intersection point
      TVector3 int_tv(intersections.at(0).ToTLorentzVector().Vect());
      length += ( box.Contain(p0) ? (v0 - int_tv).Mag() : (v1 - int_tv).Mag() ); 
    }
    // none contained -- either must have zero or two intersections
    if (n_contained == 0) {
      auto intersections = algo.Intersection(line, box);
      assert(intersections.size() == 0 || intersections.size() == 2);
      if (intersections.size() == 2) {
        TVector3 start(intersections.at(0).ToTLorentzVector().Vect());
        TVector3 end(intersections.at(1).ToTLorentzVector().Vect());
        length += (start - end).Mag();
      }
    }
  }

  return length;
}

double visibleEnergy(const simb::MCTruth &mctruth, const std::vector<sim::MCTrack> &mctrack_list, const std::vector<sim::MCShower> &mcshower_list, 
    double track_threshold, double shower_threshold) {
  double visible_E = 0;

  // total up visible energy from tracks...
  for (auto const &mct: mctrack_list) {
    // ignore neutral particles
    if (isFromNuVertex(mctruth, mct) && abs(PDGCharge(mct.PdgCode())) > 1e-4) {
      double mass = (mct.PdgCode() == 13 || mct.PdgCode() == 11) ? 0:PDGMass(mct.PdgCode());
      double this_visible_energy = mct.Start().E() - mass;
      if (this_visible_energy > track_threshold) {
        visible_E += this_visible_energy;
      }
    }
  }

  // ...and showers
  for (auto const &mcs: mcshower_list) {
    if (isFromNuVertex(mctruth, mcs) && abs(PDGCharge(mcs.PdgCode())) > 1e-4) {
      double mass = (mcs.PdgCode() == 13 || mcs.PdgCode() == 11) ? 0:PDGMass(mcs.PdgCode());
      double this_visible_energy = mcs.Start().E() - mass;
      if (this_visible_energy > shower_threshold) {
        visible_E += this_visible_energy;
      }
    }
  }

  return visible_E;
}

// define global static const PDGTable to be used by helper functions
static const TDatabasePDG *PDGTable(new TDatabasePDG);

double PDGMass(int pdg) {
  // regular particle
  if (pdg < 1000000000) {
    TParticlePDG* ple = PDGTable->GetParticle(pdg);
    return ple->Mass() * 1000.0;
  }
  // ion
  else {
    int p = (pdg % 10000000) / 10000;
    int n = (pdg % 10000) / 10 - p;
    return (PDGTable->GetParticle(2212)->Mass() * p +
            PDGTable->GetParticle(2112)->Mass() * n) * 1000.0;
  }
}

double PDGCharge(int pdg) {
  // regular particle
  if (pdg < 1000000000) {
    TParticlePDG* ple = PDGTable->GetParticle(pdg);
    return ple->Charge();
  }
  // ion
  else {
    int p = (pdg % 10000000) / 10000;
    return p * 3;
  }
}

bool isFromNuVertex(const simb::MCTruth& mc, const sim::MCShower& show, float distance)  {
  TLorentzVector nuVtx = mc.GetNeutrino().Nu().Trajectory().Position(0);
  TLorentzVector showStart = show.Start().Position();
  return (showStart - nuVtx).Mag() < distance;
}

bool isFromNuVertex(const simb::MCTruth& mc, const sim::MCTrack& track, float distance) {
  TLorentzVector nuVtx = mc.GetNeutrino().Nu().Trajectory().Position(0);
  TLorentzVector trkStart = track.Start().Position();
  return (trkStart - nuVtx).Mag() < distance;
}

  }  // namespace SBNOsc
}  // namespace ana

