#include "Constants.h"

#include <cassert>
#include <iostream>

namespace evgen {
namespace ldm {

// Useful computation
double twobody_momentum(double parent_mass, double childA_mass, double childB_mass) {
  if (parent_mass < childA_mass + childB_mass) return -1.;

  return sqrt(parent_mass * parent_mass * parent_mass * parent_mass 
    -2 * parent_mass * parent_mass * childA_mass * childA_mass
    -2 * parent_mass * parent_mass * childB_mass * childB_mass
       + childA_mass * childA_mass * childA_mass * childA_mass 
       + childB_mass * childB_mass * childB_mass * childB_mass 
    -2 * childA_mass * childA_mass * childB_mass * childB_mass) / ( 2 * parent_mass );

}

// (Somewhat) Adapted from the EnuWgt code in dk2nugenie for neutrinos.
//
// Changed to allow for a massive particle, as well as cleaned up
//
// Given the momentum in the rest frame of the parent particle, the 
// desired direction in the lab frame, and the boost vector 
// of the parent particle in the lab frame: computes the momentum and "weight"
// of the daughter going through the input end point.
//
// The weight is given as the ratio of the solid angle between the lab frame and 
// parent rest frame. I.e. -- it should be converted to an event weight as:
//
// event weight = weight * (detector lab frame solid angle) / (4*pi)
//
// For an isotropic decay in the parent rest frame
//
// Returns 0 if the calculation works. Returns -1 if the perscribed kinematics are not possible.
int calcPrtlRayWgt(double rest_frame_p, double M, TVector3 dir, TVector3 boost, double rand,
                     double& lab_frame_p_out, double& costh_rest_out, double& wgt)
{
  // preset values
  lab_frame_p_out = 0.;
  costh_rest_out = 0.;
  wgt = 0.;

  double beta = boost.Mag();
  double gamma = 1. / sqrt(1 - beta*beta);
  double E_rest = sqrt(rest_frame_p*rest_frame_p + M*M);

  // cos angle between kaon direction and daughter direction (in lab frame)
  double costh_lab = boost.Unit().Dot(dir.Unit());
  double sinth_lab_sq = 1 - costh_lab*costh_lab;

  // The problem: we are given the rest frame momentum and lab frame angle. Given this,
  // determine the lab frame momentum and rest frame angle.
  //
  // Unlike in the massless case, there are 1 OR 0 OR 2 possible solutions to this problem.
  auto costh_rest = [costh_lab, sinth_lab_sq, beta, gamma, M, rest_frame_p, E_rest](int SIGN) { return -(-SIGN * costh_lab*\
                  sqrt(-M*M + (E_rest*E_rest)/(gamma*gamma) + (M*M)*(beta*beta)*(costh_lab*costh_lab)) \
                 +E_rest*beta*gamma*sinth_lab_sq) / \
                 (rest_frame_p*gamma*(1 - beta*beta*costh_lab*costh_lab)); };

  double costh_rest_plus = costh_rest(1);
  double costh_rest_minus = costh_rest(-1);

  auto lab_frame_p = [costh_lab, sinth_lab_sq, beta, gamma, M, rest_frame_p, E_rest](int SIGN) { return (1. / (1 - beta*beta * costh_lab*costh_lab)) *\
        (E_rest * beta * costh_lab / gamma \
        + SIGN * sqrt(-M*M + (E_rest*E_rest)/(gamma*gamma) + M*M * beta*beta * costh_lab*costh_lab)); };

  double lab_frame_p_plus = lab_frame_p(1);
  double lab_frame_p_minus = lab_frame_p(-1);

  bool plus_valid = !isnan(lab_frame_p_plus) && lab_frame_p_plus > 0.;
  bool minus_valid = !isnan(lab_frame_p_minus) && lab_frame_p_minus > 0.;

  // Now compute the weight
  auto wgtf = [costh_lab, rest_frame_p, E_rest, M, beta, gamma](double p_lab, int SIGN) {
    double E_lab = sqrt(p_lab*p_lab + M*M);

    double dp_labdp_rest = (rest_frame_p / (E_rest * gamma)) * (beta*costh_lab +\
                SIGN*E_rest / (gamma * sqrt(E_rest*E_rest / (gamma*gamma) - M*M + M*M * beta*beta * costh_lab*costh_lab))) /\
                (1 - beta*beta * costh_lab*costh_lab);
    return (E_rest / E_lab) * (p_lab*p_lab / (rest_frame_p*rest_frame_p)) * abs(dp_labdp_rest);
  };

  double plus_wgt = plus_valid ? wgtf(lab_frame_p_plus, 1) : 0.;
  double minus_wgt = minus_valid ? wgtf(lab_frame_p_minus, -1) : 0.;

  //double threshold = (plus_valid || minus_valid) ? minus_valid / (plus_valid + minus_valid) : -1; 
  double threshold = 0.5;

  // Select the solution to use
  bool select_plus = plus_valid && (!minus_valid || (rand >= threshold));
  bool select_minus = minus_valid && (!plus_valid || (rand < threshold));

  // Prints out stuff
  std::cout << "COSTH LAB: " << costh_lab << std::endl;
  std::cout << "PREST: " << rest_frame_p << std::endl;
  std::cout << "MASS: " << M << std::endl;
  std::cout << "BETA: " << beta << std::endl;
  std::cout << "GAMMA: " << gamma << std::endl;
  
  std::cout << "COSTH REST PLUS: " << costh_rest_plus << std::endl;
  std::cout << "COSTH REST MINUS: " << costh_rest_minus << std::endl;
  std::cout << "PLAB PLUS: " << lab_frame_p_plus << std::endl;
  std::cout << "PLAB MINUS: " << lab_frame_p_minus << std::endl;
  std::cout << "WGT PLUS: " << plus_wgt << std::endl;
  std::cout << "WGT MINUS: " << minus_wgt << std::endl;

  if (select_plus) std::cout << "SELECTED PLUS" << std::endl;
  else if (select_minus) std::cout << "SELECTED MINUS" << std::endl;
  else std::cout << "SELECTED NONE" << std::endl;

  if (select_plus) {
    costh_rest_out = costh_rest_plus;
    lab_frame_p_out = lab_frame_p_plus;
    wgt = plus_wgt * (plus_valid + minus_valid);
  }
  else if (select_minus) {
    costh_rest_out = costh_rest_minus;
    lab_frame_p_out = lab_frame_p_minus;
    wgt = minus_wgt * (plus_valid + minus_valid);
  }
  else { // No solution!
    return -1;
  }

  return 0;
}

double minKinematicCosTheta(double parentM, double secM, double prtlM, double parentE) {
  // guard against parent-at-rest
  if (parentE - parentM < 1e-4) return 0.;

  // look up the rest frame prtl momentum
  double rest_prtlP = twobody_momentum(parentM, secM, prtlM);

  // set the dir and boost in the same direction
  double gamma = parentE / parentM;
  double beta = sqrt(1. - 1. / (gamma*gamma));

  double rest_prtlE = sqrt(prtlM*prtlM + rest_prtlP*rest_prtlP);

  return sqrt(prtlM*prtlM - rest_prtlE * rest_prtlE / (gamma*gamma)) / (prtlM*beta);
}

// Get the energy of a forward going prtl from a parent with mass M / energy E
double forwardPrtlEnergy(double parentM, double secM, double prtlM, double parentE) {
  // look up the rest frame prtl momentum
  double rest_prtlP = twobody_momentum(parentM, secM, prtlM);

  // set the dir and boost in the same direction
  double gamma = parentE / parentM;
  double beta = sqrt(1. - 1. / (gamma*gamma));
  TVector3 dir(0., 0., 1.);
  TVector3 boost(0., 0., beta);

  // make sure we get the positive (more energetic) solution
  double rand = 1.;

  double lab_frame_p_out, costh_rest_out, wgt;
  int err = calcPrtlRayWgt(rest_prtlP, prtlM, dir, boost, rand,
    lab_frame_p_out, costh_rest_out, wgt);
  assert(!err);
  (void) err;

  double lab_prtlE = sqrt(lab_frame_p_out*lab_frame_p_out + prtlM*prtlM);
  return lab_prtlE;
}

// TODO: use particle database
double secPDG2Mass(int pdg) {
  switch (pdg) {
    case 11:
    case -11:
      return elec_mass;
      break;
    case 13:
    case -13:
      return muon_mass;
      break;
    case 211:
    case -211:
      return pionp_mass;
      break;
    default:
      std::cerr << "BAD secondary pdg: " << pdg << std::endl;
      std::cout << "BAD secondary pdg: " << pdg << std::endl;
      return 0.;
  }
}

} // end namespace
} // end namespace