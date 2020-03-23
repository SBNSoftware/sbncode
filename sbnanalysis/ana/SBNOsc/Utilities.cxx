#include <iostream>

#include "TDatabasePDG.h"

#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCNeutrino.h"
#include "core/Event.hh"

#include "Utilities.h"

#include "ubcore/LLBasicTool/GeoAlgo/GeoAABox.h"
#include "ubcore/LLBasicTool/GeoAlgo/GeoAlgo.h"
#include "ubcore/LLBasicTool/GeoAlgo/GeoLineSegment.h"

#include <TMath.h>
#include <TRandom3.h>

namespace ana {
  namespace SBNOsc {

void hello() {
  std::cout << "Hello SBNOsc!" << std::endl;
}


event::Interaction TruthReco(const simb::MCTruth& mctruth) {
  event::Interaction interaction;

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
    event::FinalStateParticle fsp;
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


double ECCQE(const TVector3& l_momentum, double l_energy,
             double energy_distortion, double angle_distortion) {
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


double NuMuOscillation(double numu_energy, double numu_dist,
                       double osc_dm2, double osc_angle) {
  double overlap = sin(2*osc_angle);
  double energy_factor = sin(1.27 * osc_dm2 * numu_dist / numu_energy);
  return 1 - overlap * overlap * energy_factor * energy_factor;
}


double containedLength(const TVector3 &v0, const TVector3 &v1,
                       const std::vector<geoalgo::AABox> &boxes) {
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
      // Because of floating point errors, it can sometimes happen
      // that there is 1 contained point but no "Intersections"
      // if one of the points is right on the edge
      if (intersections.size() == 0) {
        // determine which point is on the edge
        double tol = 1e-5;
        bool p0_edge = algo.SqDist(p0, box) < tol;
        bool p1_edge = algo.SqDist(p1, box) < tol;
        assert(p0_edge || p1_edge);
        // contained one is on edge -- can treat both as not contained
        //
        // In this case, no length
        if ((p0_edge && box.Contain(p0)) || (box.Contain(p1) && p1_edge))
          continue;
        // un-contaned one is on edge -- treat both as contained
        else if ((p0_edge && box.Contain(p1)) || (box.Contain(p0) && p1_edge)) {
	  length = (v1 - v0).Mag();
	  break;
        }
        else {
          assert(false); // bad
        }
      }
      // floating point errors can also falsely cause 2 intersection points
      //
      // in this case, one of the intersections must be very close to the 
      // "contained" point, so the total contained length will be about
      // the same as the distance between the two intersection points
      else if (intersections.size() == 2) {
        length += (intersections.at(0).ToTLorentzVector().Vect() - intersections.at(1).ToTLorentzVector().Vect()).Mag();
        continue;
      }
      // "Correct"/ideal case -- 1 intersection point
      else if (intersections.size() == 1) {
        // get TVector at intersection point
        TVector3 int_tv(intersections.at(0).ToTLorentzVector().Vect());
        length += ( box.Contain(p0) ? (v0 - int_tv).Mag() : (v1 - int_tv).Mag() ); 
      }
      else assert(false); // bad
    }
    // none contained -- either must have zero or two intersections
    if (n_contained == 0) {
      auto intersections = algo.Intersection(line, box);
      if (!(intersections.size() == 0 || intersections.size() == 2)) {
        // more floating point error fixes...
        //
        // figure out which points are near the edge
        double tol = 1e-5;
        bool p0_edge = algo.SqDist(p0, box) < tol;
        bool p1_edge = algo.SqDist(p1, box) < tol;
        // and which points are near the intersection
        TVector3 vint = intersections.at(0).ToTLorentzVector().Vect();
        bool p0_int = (v0 - vint).Mag() < tol;
        bool p1_int = (v1 - vint).Mag() < tol;
        // exactly one of them should produce the intersection
        assert((p0_int && p0_edge) != (p1_int && p1_edge));
        // both close to edge -- full length is contained
        if (p0_edge && p1_edge) {
          length += (v0 - v1).Mag();
        }
        // otherwise -- one of them is not on an edge, no length is contained
        else {}
      }
      // assert(intersections.size() == 0 || intersections.size() == 2);
      else if (intersections.size() == 2) {
        TVector3 start(intersections.at(0).ToTLorentzVector().Vect());
        TVector3 end(intersections.at(1).ToTLorentzVector().Vect());
        length += (start - end).Mag();
      }
    }
  }

  return length;
}

// run the proposal energy smearing w/ MCParticles
double visibleEnergyProposalMCParticles(TRandom &rand, const simb::MCTruth &mctruth, const std::vector<sim::MCTrack> mctrack_list, const VisibleEnergyCalculator &calculator) {
  double total = 0;
  for (int iparticle=0; iparticle<mctruth.NParticles(); iparticle++) {
    const simb::MCParticle& particle = mctruth.GetParticle(iparticle);

    int pdg = particle.PdgCode();

    // if there's no mass entry, we probably don't care about it
    if (PDGMass(pdg) < 0) continue;

    double track_energy = (particle.Momentum().E() - PDGMass(pdg) / 1000. /* MeV -> GeV */);
    double smear_energy = track_energy;

    // In the proposal, different PDG codes got different energies

    // kaons did not have the mass subtracted
    if (abs(pdg) == 321) {
      track_energy = particle.Momentum().E();
      smear_energy = track_energy;
    }
    // pions did not have mass subtracted in the smearing
    if (abs(pdg) == 211) {
      smear_energy = particle.Momentum().E();
    }

    bool skip = false;

    // threshold -- only for protons
    if (abs(pdg) == 2212 && track_energy < calculator.track_threshold) {
      skip = true;
    } 

    // add in smeared energy
    double this_smeared_energy = rand.Gaus(track_energy, smear_energy * calculator.track_energy_distortion);
    // clamp to zero
    this_smeared_energy = std::max(this_smeared_energy, 0.);

    // Ignore particles not from nu vertex
    if (!isFromNuVertex(mctruth, particle) || particle.Process() != "primary") {
      skip = true;
    }
    // ignore everything except for protons and kaons and pions
    if (!(abs(pdg) == 2212 || abs(pdg) == 211 || abs(pdg) == 321)) {
      skip = true;
    }
    // double distance = (mctruth.GetNeutrino().Nu().Trajectory().Position(0).Vect() - mct.Start().Position().Vect()).Mag();
    // double l_distance = (mctruth.GetNeutrino().Nu().Trajectory().Position(0) - mct.Start().Position()).Mag();
    // std::cout << "New particle -- true E: " << track_energy << " smeared E: " << this_smeared_energy << " PDG: " << pdg << " is primary: " << (calculator.lepton_index == ind) << " is skipped: " << skip << " distance: " << distance << " lorentz distance: " << l_distance << " process: " << mct.Process() << std::endl; 
    if (skip) continue;

    total += this_smeared_energy;
  }

  // ...and primary lepton energy (for CC events)
  // only add in extra here if identified "lepton" is actually a lepton
  if (calculator.lepton_index >= 0 && (abs(mctrack_list[calculator.lepton_index].PdgCode()) == 13 || abs(mctrack_list[calculator.lepton_index].PdgCode()) == 11)) {
    double lepton_energy = smearLeptonEnergy(rand, mctrack_list[calculator.lepton_index], calculator);
    total += lepton_energy;
  }
  // std::cout << "Reco Energy: " << total << std::endl;

  return total;
  
  }

// copies the exact energy smearing used in the proposal
double visibleEnergyProposal(TRandom &rand, const simb::MCTruth &mctruth, const std::vector<sim::MCTrack> &mctrack_list, const VisibleEnergyCalculator &calculator) {
  double total = 0.;

  // std::cout << "\n\nNew Interaction" << std::endl;
  // std::cout << "True E: " << mctruth.GetNeutrino().Nu().E() << std::endl;
  for (unsigned ind = 0; ind < mctrack_list.size(); ind++) {
    auto const &mct = mctrack_list[ind];
    int pdg = mct.PdgCode();

    double track_energy = (mct.Start().E() - PDGMass(pdg)) / 1000. /* MeV -> GeV */;
    double smear_energy = track_energy;

    // In the proposal, different PDG codes got different energies

    // kaons did not have the mass subtracted
    if (abs(pdg) == 321) {
      track_energy = mct.Start().E() / 1000.;
      smear_energy = track_energy;
    }
    // pions did not have mass subtracted in the smearing
    if (abs(pdg) == 211) {
      smear_energy = mct.Start().E() / 1000.;
    }

    bool skip = false;

    // threshold -- only for protons
    if (abs(pdg) == 2212 && track_energy < calculator.track_threshold) {
      skip = true;
    } 

    // add in smeared energy
    double this_smeared_energy = rand.Gaus(track_energy, smear_energy * calculator.track_energy_distortion);
    // clamp to zero
    this_smeared_energy = std::max(this_smeared_energy, 0.);

    // Ignore particles not from nu vertex
    if (!isFromNuVertex(mctruth, mct) || mct.Process() != "primary") {
      skip = true;
    }
    // ignore everything except for protons and kaons and pions
    if (!(abs(pdg) == 2212 || abs(pdg) == 211 || abs(pdg) == 321)) {
      skip = true;
    }
    // account for primary track later
    if ((abs(pdg) == 13 || abs(pdg) == 11) && calculator.lepton_index == ind) {
      skip = true;
    }
    // double distance = (mctruth.GetNeutrino().Nu().Trajectory().Position(0).Vect() - mct.Start().Position().Vect()).Mag();
    // double l_distance = (mctruth.GetNeutrino().Nu().Trajectory().Position(0) - mct.Start().Position()).Mag();
    // std::cout << "New particle -- true E: " << track_energy << " smeared E: " << this_smeared_energy << " PDG: " << pdg << " is primary: " << (calculator.lepton_index == ind) << " is skipped: " << skip << " distance: " << distance << " lorentz distance: " << l_distance << " process: " << mct.Process() << std::endl; 
    if (skip) continue;

    total += this_smeared_energy;
  }

  // ...and primary lepton energy (for CC events)
  // only add in extra here if identified "lepton" is actually a lepton
  if (calculator.lepton_index >= 0 && (abs(mctrack_list[calculator.lepton_index].PdgCode()) == 13 || abs(mctrack_list[calculator.lepton_index].PdgCode()) == 11)) {
    double lepton_energy = smearLeptonEnergy(rand, mctrack_list[calculator.lepton_index], calculator);
    total += lepton_energy;
  }
  // std::cout << "Reco Energy: " << total << std::endl;

  return total;
}

std::pair<double, double> visibleEnergySplit(TRandom &rand, 
    const simb::MCTruth &mctruth, 
    const std::vector<sim::MCTrack> &mctrack_list, 
    const VisibleEnergyCalculator &calculator) {

  double visible_E = 0;

  // primary leptron track
  const sim::MCTrack *lepton_track = NULL;
  bool lepton_track_exists = false;

  // total up visible energy from tracks...
  double track_visible_energy = 0.;
  for (unsigned ind = 0; ind < mctrack_list.size(); ind++) {
    auto const &mct = mctrack_list[ind];
    // ignore particles not from nu vertex, non primary particles, and uncharged particles
    if (!isFromNuVertex(mctruth, mct) || abs(PDGCharge(mct.PdgCode())) < 1e-4 || mct.Process() != "primary")
       continue;
    // account for primary lepton later
    if ((abs(mct.PdgCode()) == 13 || abs(mct.PdgCode()) == 11) && calculator.lepton_index == ind) {
      continue;
    }
    // account for primary pion later
    if ((abs(mct.PdgCode()) == 211) && calculator.lepton_index == ind) {
      continue;
    }

    double mass = PDGMass(mct.PdgCode());
    double this_visible_energy = (mct.Start().E() - mass) / 1000. /* MeV to GeV */;
    if (this_visible_energy > calculator.track_threshold) {
      track_visible_energy += this_visible_energy;
    }
  }

  // do energy smearing
  if (calculator.track_energy_distortion > 1e-4) {
    track_visible_energy = rand.Gaus(track_visible_energy, track_visible_energy*calculator.track_energy_distortion);
    // clamp to 0
    track_visible_energy = std::max(track_visible_energy, 0.);
  }

  double hadronic_visible_E = track_visible_energy;

  double leptonic_visible_E = 0.;

  // ...and primary lepton energy (for CC events)
  // only add in extra here if identified "lepton" is actually a lepton
  if (calculator.lepton_index >= 0 && (abs(mctrack_list[calculator.lepton_index].PdgCode()) == 13 || abs(mctrack_list[calculator.lepton_index].PdgCode()) == 11)) {
    leptonic_visible_E = smearLeptonEnergy(rand, mctrack_list[calculator.lepton_index], calculator);
  }
  else if (calculator.lepton_index >= 0 && (abs(mctrack_list[calculator.lepton_index].PdgCode()) == 211)) {
    const sim::MCTrack &mct = mctrack_list[calculator.lepton_index];
    double pion_visible_energy =  (mct.Start().E() -PDGMass(mct.PdgCode())) / 1000. /* MeV to GeV */;
    leptonic_visible_E = std::max(rand.Gaus(track_visible_energy, track_visible_energy*calculator.track_energy_distortion), 0.);
  }

  return {hadronic_visible_E, leptonic_visible_E};
}

double visibleEnergy(TRandom &rand, const simb::MCTruth &mctruth, const std::vector<sim::MCTrack> &mctrack_list, const std::vector<sim::MCShower> &mcshower_list, 
    const VisibleEnergyCalculator &calculator, bool include_showers) {
  double visible_E = 0;

  // primary leptron track
  const sim::MCTrack *lepton_track = NULL;
  bool lepton_track_exists = false;

  // total up visible energy from tracks...
  double track_visible_energy = 0.;
  for (unsigned ind = 0; ind < mctrack_list.size(); ind++) {
    auto const &mct = mctrack_list[ind];
    // ignore particles not from nu vertex, non primary particles, and uncharged particles
    if (!isFromNuVertex(mctruth, mct) || abs(PDGCharge(mct.PdgCode())) < 1e-4 || mct.Process() != "primary")
       continue;
    // account for primary lepton later
    if ((abs(mct.PdgCode()) == 13 || abs(mct.PdgCode()) == 11) && calculator.lepton_index == ind) {
      continue;
    }
    // ignore non-primary pion??
    //if ((abs(mct.PdgCode()) == 211) && calculator.lepton_index != ind) {
    //  continue;
    //}

    double mass = PDGMass(mct.PdgCode());
    double this_visible_energy = (mct.Start().E() - mass) / 1000. /* MeV to GeV */;
    if (this_visible_energy > calculator.track_threshold) {
      track_visible_energy += this_visible_energy;
    }
  }

  // do energy smearing
  if (calculator.track_energy_distortion > 1e-4) {
    track_visible_energy = rand.Gaus(track_visible_energy, track_visible_energy*calculator.track_energy_distortion);
    // clamp to 0
    track_visible_energy = std::max(track_visible_energy, 0.);
  }

  // ...and showers
  double shower_visible_energy = 0.;
  if (include_showers) {
    for (auto const &mcs: mcshower_list) {
      // ignore particles not from nu vertex, non primary particles, and uncharged particles
      if (!isFromNuVertex(mctruth, mcs) || abs(PDGCharge(mcs.PdgCode())) < 1e-4 || mcs.Process() != "primary")
        continue; 
      // account for primary lepton later
      if ((abs(mcs.PdgCode()) == 13 || abs(mcs.PdgCode()) == 11) && isFromNuVertex(mctruth, mcs))
        continue; 

      double mass = PDGMass(mcs.PdgCode());
      double this_visible_energy = (mcs.Start().E() - mass) / 1000. /* MeV to GeV */;
      
      if (this_visible_energy > calculator.shower_threshold) {
        shower_visible_energy += this_visible_energy;
      }
    }
  }

  // do energy smearing
  if (calculator.shower_energy_distortion > 1e-4) {
    shower_visible_energy = rand.Gaus(shower_visible_energy, shower_visible_energy*calculator.shower_energy_distortion);
    // clamp to 0
    shower_visible_energy = std::max(shower_visible_energy, 0.);
  }

  visible_E = track_visible_energy + shower_visible_energy;

  // ...and primary lepton energy (for CC events)
  // only add in extra here if identified "lepton" is actually a lepton
  if (calculator.lepton_index >= 0 && (abs(mctrack_list[calculator.lepton_index].PdgCode()) == 13 || abs(mctrack_list[calculator.lepton_index].PdgCode()) == 11)) {
    visible_E += smearLeptonEnergy(rand, mctrack_list[calculator.lepton_index], calculator);
  }

  return visible_E;
}


double smearLeptonEnergy(TRandom &rand, const sim::MCTrack &mct, const VisibleEnergyCalculator &calculator) {
  // if not contained and zero length, return
  if (!calculator.lepton_contained && calculator.lepton_contained_length < 1e-4) {
    // std::cout << "Out of detector lepton\n";
    return 0;
  }

  double smearing_percentage;
  if (calculator.lepton_contained) {
    smearing_percentage = calculator.lepton_energy_distortion_contained;
  } else {
    double A = calculator.lepton_energy_distortion_leaving_A;
    double B = calculator.lepton_energy_distortion_leaving_B;
    smearing_percentage = -A * TMath::Log(B * calculator.lepton_contained_length);
  }
  // smear visible energy
  double lepton_visible_energy = (mct.Start().E()) / 1000.; /* MeV to GeV */
  double smeared_lepton_visible_energy = rand.Gaus(lepton_visible_energy, smearing_percentage * lepton_visible_energy);
  // clamp to 0
  smeared_lepton_visible_energy = std::max(smeared_lepton_visible_energy, 0.);

  // std::cout << "Lepton -- is_contained: " << calculator.lepton_contained << " length: " << calculator.lepton_contained_length << " smearing: " << smearing_percentage << " true E: " << lepton_visible_energy << " smeared E: " << smeared_lepton_visible_energy << std::endl; 

  return smeared_lepton_visible_energy;
}


// define global static const PDGTable to be used by helper functions
static const TDatabasePDG *PDGTable(new TDatabasePDG);


double PDGMass(int pdg) {
  // regular particle
  if (pdg < 1000000000) {
    TParticlePDG* ple = PDGTable->GetParticle(pdg);
    if (ple == NULL) return -1;
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

bool isFromNuVertex(const simb::MCTruth& mc, const simb::MCParticle& mcp, float distance) {
  TVector3 nuVtx = mc.GetNeutrino().Nu().Position().Vect();
  TVector3 partStart = mcp.Position().Vect();
  return (partStart - nuVtx).Mag() < distance;
}

bool isFromNuVertex(const simb::MCTruth& mc, const sim::MCShower& show,
                    float distance)  {
  TVector3 nuVtx = mc.GetNeutrino().Nu().Trajectory().Position(0).Vect();
  TVector3 showStart = show.Start().Position().Vect();
  return (showStart - nuVtx).Mag() < distance;
}


bool isFromNuVertex(const simb::MCTruth& mc, const sim::MCTrack& track,
                    float distance) {
  TVector3 nuVtx = mc.GetNeutrino().Nu().Trajectory().Position(0).Vect();
  TVector3 trkStart = track.Start().Position().Vect();
  return (trkStart - nuVtx).Mag() < distance;
}

  }  // namespace SBNOsc
}  // namespace ana
