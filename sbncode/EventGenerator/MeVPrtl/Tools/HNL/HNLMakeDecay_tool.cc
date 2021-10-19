/**
 *
 */

// Framework Includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Utilities/ToolMacros.h"
#include "cetlib/cpu_timer.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "sbnobj/Common/EventGen/MeVPrtl/MeVPrtlFlux.h"

// local includes
#include "sbncode/EventGenerator/MeVPrtl/Tools/IMeVPrtlDecay.h"
#include "sbncode/EventGenerator/MeVPrtl/Tools/Constants.h"

// LArSoft includes
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"
#include "nusimdata/SimulationBase/MCParticle.h"

// std includes
#include <string>
#include <iostream>
#include <memory>
#include <utility>

// math
#include <math.h>
#include <gsl/gsl_integration.h>

// constants
#include "TDatabasePDG.h"

#include "HNLDecayDalitz.h"

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace evgen {
namespace ldm {
/**
 *  @brief  HNLMakeDecay class definiton
 *
 *  Implementation of HNL decay ->mupi taken from:
 *      https://arxiv.org/abs/1610.08512
 */
class HNLMakeDecay : public IMeVPrtlDecay {
public:
    /**
     *  @brief  Constructor
     */
    HNLMakeDecay(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~HNLMakeDecay();

    void configure(fhicl::ParameterSet const &pset) override;

    bool Decay(const MeVPrtlFlux &flux, const TVector3 &in, const TVector3 &out, MeVPrtlDecay &decay, double &weight) override;

    // returns the max weight of configured
    float MaxWeight() override { 
      return fMaxWeight; 
    }

private:
  // Internal data
  double fMaxWeight;
  gsl_integration_workspace *fIntegrator;
  unsigned fIntegratorSize;

  // Configure the MaxWeight
  double fReferenceUE4;
  double fReferenceUM4;
  double fReferenceHNLMass;
  double fReferenceRayLength;
  double fReferenceHNLEnergy;
  double fReferenceHNLKaonEnergy;

  // Guardrail for small decay lengths
  double fMinDetectorDistance;

  // Configure the particle
  bool fMajorana;

  // Internal struct for holding decay information
  struct DecayFinalState {
    float width;
    std::vector<TLorentzVector> mom;
    std::vector<int> pdg;
  };

  // In the threebody-decay case, we need to specify the three momentum vectors, not just the overall
  // magnitude
  struct ThreebodyMomentum {
    TLorentzVector A;
    TLorentzVector B;
    TLorentzVector C;
  };

  typedef DecayFinalState(HNLMakeDecay::*HNLDecayFunction)(const MeVPrtlFlux &flux);

  std::map<std::string, HNLDecayFunction> fAvailableDecays;
  std::map<std::string, float> fAvailableDecayMasses;
  std::vector<std::string> fDecayConfig;
  std::vector<HNLDecayFunction> fSelectedDecays;

  // Helper functions
  double CalculateMaxWeight();
  double CalculateKDARDecayLength();
  ThreebodyMomentum isotropic_threebody_momentum(double parent_mass, double childA_mass, double childB_mass, double childC_mass);
  double I1(double x, double y, double z);
  double I2(double x, double y, double z);
  double NuDiLepDecayWidth(double hnl_mass, double ue4, double um4, bool nu_is_muon, bool lep_is_muon);

  // Decay implementation functions
  DecayFinalState NuMupMum(const MeVPrtlFlux &flux);
  DecayFinalState LepPi(const MeVPrtlFlux &flux, bool is_muon);
  DecayFinalState EPi(const MeVPrtlFlux &flux) { return LepPi(flux, false); }
  DecayFinalState MuPi(const MeVPrtlFlux &flux) { return LepPi(flux, true); }
};

// helpers
double lambda(double a, double b, double c) {
  return a*a + b*b + c*c - 2*a*b - 2*b*c - 2*c*a;
}


// Valid for decays where the matix element has no kinematic dependence (i.e. a constant Dalitz density)
HNLMakeDecay::ThreebodyMomentum HNLMakeDecay::isotropic_threebody_momentum(double parent_mass, double childA_mass, double childB_mass, double childC_mass) {
  ThreebodyMomentum ret;
  double sumofdaughtermass = childA_mass + childB_mass + childC_mass;
  if (parent_mass < sumofdaughtermass) { // shouldn't happen
    return ret;
  } 

  double E_A, E_B, E_C;
  double P_A, P_B, P_C;
  double P_max, P_sum;
  // Kinetic Energy is distributed uniformly to daughters, with the constraint that
  // total momentum is conserved. Randomly allocate energy until a possible such configuration
  // is found
  do {
    double r1 = GetRandom();
    double r2 = GetRandom();
  
    E_A = childA_mass + (parent_mass - sumofdaughtermass) * std::min(r1, r2);
    E_B = childB_mass + (parent_mass - sumofdaughtermass) * std::min(1-r1,1-r2);
    E_C = childC_mass + (parent_mass - sumofdaughtermass) * abs(r1-r2);

    P_A = sqrt(E_A*E_A - childA_mass*childA_mass);
    P_B = sqrt(E_B*E_B - childB_mass*childB_mass);
    P_C = sqrt(E_C*E_C - childC_mass*childC_mass);

    P_max = std::max(std::max(P_A,P_B),P_C);
    P_sum = P_A + P_B + P_C;

  } while(P_max > P_sum - P_max);

  // Found a valid momentum allocation!
  
  // Pick a random direction for A, have the direction of B, C work to conserve momentum
  TVector3 dirA = RandomUnitVector();

  // daughter particles B and C have the same momentum perpindicular to the direction of A
  // Solving for the direction along the axis of particle A gives:
  double cos_thAB = (P_C*P_C - P_B*P_B - P_A*P_A) / (2. * P_A * P_B);
  double sin_thAB = sqrt(1. - cos_thAB * cos_thAB);
  double cos_thAC = (P_B*P_B - P_C*P_C - P_A*P_A) / (2. * P_A * P_C);
  double sin_thAC = sqrt(1. - cos_thAC * cos_thAC);

  // The azimuthal angle of B and C about A is distributed uniformly
  double gammaB = (2*GetRandom() - 1.) * M_PI;
  double gammaC = -gammaB;

  TVector3 dirB(
  sin_thAB*cos(gammaB)*dirA.CosTheta()*sin(dirA.Phi()) - sin_thAB*sin(gammaB)*sin(dirA.Phi()) + cos_thAB*sqrt(1.-dirA.CosTheta()*dirA.CosTheta()) * cos(dirA.Phi()),
  sin_thAB*cos(gammaB)*dirA.CosTheta()*cos(dirA.Phi()) - sin_thAB*sin(gammaB)*cos(dirA.Phi()) + cos_thAB*sqrt(1.-dirA.CosTheta()*dirA.CosTheta()) * sin(dirA.Phi()),
 -sin_thAB*cos(gammaB)*sqrt(1. - dirA.CosTheta() * dirA.CosTheta()) + cos_thAB*dirA.CosTheta());

  TVector3 dirC(
  sin_thAC*cos(gammaC)*dirA.CosTheta()*sin(dirA.Phi()) - sin_thAC*sin(gammaC)*sin(dirA.Phi()) + cos_thAC*sqrt(1.-dirA.CosTheta()*dirA.CosTheta()) * cos(dirA.Phi()),
  sin_thAC*cos(gammaC)*dirA.CosTheta()*cos(dirA.Phi()) - sin_thAC*sin(gammaC)*cos(dirA.Phi()) + cos_thAC*sqrt(1.-dirA.CosTheta()*dirA.CosTheta()) * sin(dirA.Phi()),
 -sin_thAC*cos(gammaC)*sqrt(1. - dirA.CosTheta() * dirA.CosTheta()) + cos_thAC*dirA.CosTheta());

  ret.A = TLorentzVector(P_A*dirA, E_A);
  ret.B = TLorentzVector(P_B*dirB, E_B);
  ret.C = TLorentzVector(P_C*dirC, E_C);

  return ret;
}

double I1_integrand(double s, void *param) {
  double *xyz = (double *)param;
  double x = xyz[0];
  double y = xyz[1];
  double z = xyz[2];

  return 12.*(s - x*x - y*y)*(1 + z*z - s)*sqrt(lambda(s,x*x,y*y))*sqrt(lambda(1.,s,z*z))/s;
}

double I2_integrand(double s, void *param) {
  double *xyz = (double *)param;
  double x = xyz[0];
  double y = xyz[1];
  double z = xyz[2];

  return 24*y*z*(1. + x*x - s)*sqrt(lambda(s,y*y,z*z))*sqrt(lambda(1.,s,z*z))/s;
}

double HNLMakeDecay::I1(double x, double y, double z) {
  gsl_function F;
  double xyz[3];
  xyz[0] = x;
  xyz[1] = y;
  xyz[2] = z;
  F.function = &I1_integrand;
  F.params = xyz;

  double result, error;
  gsl_integration_qags(&F, (x+y)*(x+y), (1.-z)*(1.-z), 0., 1e-7, fIntegratorSize, fIntegrator, &result, &error);

  return result;
}

double HNLMakeDecay::I2(double x, double y, double z) {
  gsl_function F;
  double xyz[3];
  xyz[0] = x;
  xyz[1] = y;
  xyz[2] = z;
  F.function = &I2_integrand;
  F.params = xyz;

  double result, error;
  gsl_integration_qags(&F, (y+z)*(y+z), (1.-x)*(1.-x), 0., 1e-7, fIntegratorSize, fIntegrator, &result, &error);

  return result;
}


double HNLMakeDecay::NuDiLepDecayWidth(double hnl_mass, double ue4, double um4, bool nu_is_muon, bool lep_is_muon) {
  double hnl_mass_pow5 = hnl_mass*hnl_mass*hnl_mass*hnl_mass*hnl_mass;
  double u4 = nu_is_muon ? um4 : ue4;
  double lep_mass = lep_is_muon ? Constants::Instance().muon_mass : Constants::Instance().elec_mass;

  double Gfermi = Constants::Instance().Gfermi;
  double gL = Constants::Instance().gL;
  double gR = Constants::Instance().gR;

  if (hnl_mass < lep_mass * 2.) return 0.;

  int CC = (lep_is_muon == nu_is_muon);

  double I1val = I1(0., lep_mass / hnl_mass, lep_mass / hnl_mass);
  double I2val = I2(0., lep_mass / hnl_mass, lep_mass / hnl_mass);

  double width = (Gfermi*Gfermi*hnl_mass_pow5) * u4 * ((gL*gR/*NC*/ + CC*gR/*CC*/)*I1val + (gL*gL+gR*gR+CC*(1+2.*gL))*I2val);

  return width;
}

HNLMakeDecay::DecayFinalState HNLMakeDecay::NuMupMum(const MeVPrtlFlux &flux) {
  HNLMakeDecay::DecayFinalState ret;

  double muon_mass = Constants::Instance().muon_mass;

  // Decay not kinematically allowed
  if (2*muon_mass > flux.mass) {
    ret.width = 0.;
    return ret;
  } 

  double ue4 = flux.C1;
  double um4 = flux.C2;
  double nue_width = NuDiLepDecayWidth(flux.mass, ue4, um4, false, true);
  double numu_width = NuDiLepDecayWidth(flux.mass, ue4, um4, true, true);

  ret.width = nue_width + numu_width;

  if (ret.width == 0.) return ret;

  // Majorana gets an extra factor b.c. it can go to nu and nubar
  if (fMajorana) ret.width *= 2;

  // Three body decay
  //
  // TODO: account for anisotropies in decay
  ThreebodyMomentum momenta = isotropic_threebody_momentum(flux.mass, 0., muon_mass, muon_mass); 

  // pick whether the neutrino is nue or numu
  int nu_pdg_sign;
  if (fMajorana) {
    nu_pdg_sign = (GetRandom() > 0.5) ? 1:-1;
  }
  else {
    // same as the HNL
    nu_pdg_sign = (flux.secondary_pdg > 0) ? -1 : 1;
  }
  int nu_pdg = nu_pdg_sign * ((GetRandom() > numu_width / (numu_width + nue_width)) ? 14 : 12);
  ret.pdg.push_back(nu_pdg);
  ret.mom.push_back(momenta.A);

  ret.pdg.push_back(13);
  ret.mom.push_back(momenta.B);
  ret.pdg.push_back(-13);
  ret.mom.push_back(momenta.C);

  return ret;
}

HNLMakeDecay::DecayFinalState HNLMakeDecay::LepPi(const MeVPrtlFlux &flux, bool is_muon) {
  HNLMakeDecay::DecayFinalState ret;
  double lep_mass = is_muon ? Constants::Instance().muon_mass : Constants::Instance().elec_mass;
  int lep_pdg = is_muon ? 13 : 11;

  double piplus_mass = Constants::Instance().piplus_mass;
  double Gfermi = Constants::Instance().Gfermi;
  double fpion = Constants::Instance().fpion;
  double abs_Vud_squared = Constants::Instance().abs_Vud_squared;

  // Decay not kinematically allowed
  if (lep_mass + piplus_mass > flux.mass) {
    ret.width = 0.;
    return ret;
  }

  double u4 = is_muon ? flux.C2 : flux.C1;
  double lep_ratio = (lep_mass * lep_mass) / (flux.mass * flux.mass);
  double pion_ratio = (piplus_mass * piplus_mass) / (flux.mass * flux.mass);
  double Ifunc = ((1 + lep_ratio + pion_ratio)*(1.+lep_ratio) - 4*lep_ratio) * sqrt(lambda(1., lep_ratio, pion_ratio));
  ret.width = u4 * (Gfermi * Gfermi *fpion * fpion * abs_Vud_squared * flux.mass * flux.mass * flux.mass * Ifunc) / (16 * M_PI);

  // Majorana gets an extra factor b.c. it can go to pi+l- and pi-l+
  if (fMajorana) ret.width *= 2;

  // Majorana decays don't conserve lepton number, Dirac decay's do
  int lep_pdg_sign;
  if (fMajorana) {
    lep_pdg_sign = (GetRandom() > 0.5) ? 1 : -1;
  }
  else {
    // Dirac HNL caries opposite lepton number to production lepton
    lep_pdg_sign = (flux.secondary_pdg > 0) ? -1 : 1;
  }

  // Use rejection sampling to draw a direction for the child particles
  //
  // Work in the lab frame
  double dalitz_max = HNLLepPiDalitzMax(Constants::Instance().kplus_mass, flux.sec.M(), flux.mass, piplus_mass, lep_mass); 
  double this_dalitz = 0.;
  double p = evgen::ldm::twobody_momentum(flux.mass, lep_mass, piplus_mass);
  TLorentzVector LB;
  TLorentzVector PI;
  do {
    TVector3 dir = RandomUnitVector(); 
    LB = TLorentzVector(p*dir, sqrt(p*p + lep_mass*lep_mass));
    PI = TLorentzVector(-p*dir, sqrt(p*p + piplus_mass*piplus_mass));
    LB.Boost(flux.mom.BoostVector());
    PI.Boost(flux.mom.BoostVector());
    
    this_dalitz = ((flux.secondary_pdg > 0 ) != (lep_pdg_sign > 0)) ? \
      evgen::ldm::HNLLepPiLNCDalitz(flux.kmom, flux.sec, flux.mom, PI, LB):
      evgen::ldm::HNLLepPiLNVDalitz(flux.kmom, flux.sec, flux.mom, PI, LB);

    assert(this_dalitz < dalitz_max);
    if (this_dalitz > dalitz_max) {
      std::cerr << "VERY VERY BAD!!!! Incorrect dalitz max!!!\n";
      std::cout << "VERY VERY BAD!!!! Incorrect dalitz max!!!\n";
      std::cout << "PK: " << flux.kmom.E() << " " << flux.kmom.Px() << " " << flux.kmom.Py() << " " << flux.kmom.Pz() << std::endl;
      std::cout << "PA: " << flux.sec.E() << " " << flux.sec.Px() << " " << flux.sec.Py() << " " << flux.sec.Pz() << std::endl;
      std::cout << "PN: " << flux.mom.E() << " " << flux.mom.Px() << " " << flux.mom.Py() << " " << flux.mom.Pz() << std::endl;
      std::cout << "PP: " << PI.E() << " " << PI.Px() << " " << PI.Py() << " " << PI.Pz() << std::endl;
      std::cout << "PB: " << LB.E() << " " << LB.Px() << " " << LB.Py() << " " << LB.Pz() << std::endl;

      std::cout << "This Dalitz: " << this_dalitz << std::endl;
      std::cout << "Max Dalitz: " << dalitz_max << std::endl;
      std::cout << "LNC: " << ((flux.secondary_pdg > 0 ) != (lep_pdg_sign > 0)) << std::endl;

      exit(1);
    }
  } while (GetRandom() > this_dalitz / dalitz_max);

  // lepton
  ret.mom.push_back(LB);
  ret.pdg.push_back(lep_pdg*lep_pdg_sign);

  // pion
  ret.mom.emplace_back(PI);
  ret.pdg.push_back(211*lep_pdg_sign); // negative of lepton-charge has same-sign-PDG code

  return ret;
}

double HNLMakeDecay::CalculateKDARDecayLength() {
  double ue4 = fReferenceUE4;
  double um4 = fReferenceUM4;
  double hnl_mass = fReferenceHNLMass;
  // If the muon coupling is turned on, the lower energy / lower decay length HNL would have a muon
  // secondary
  double lep_mass = (um4 > 0.) ? Constants::Instance().muon_mass : Constants::Instance().elec_mass;
  double P = twobody_momentum(Constants::Instance().kplus_mass, lep_mass, hnl_mass);
  double E = sqrt(P*P + hnl_mass * hnl_mass);

  // make a fake HNL flux with these parameters
  MeVPrtlFlux hnl;
  hnl.C1 = ue4;
  hnl.C2 = um4;
  hnl.mass = hnl_mass;
  hnl.mom = TLorentzVector(P * TVector3(1, 0, 0), E);
  hnl.secondary_pdg = -1; // doesn't affect the end width, just make it viable

  // Mock up the secondary momenta -- this shouldn't affect width and just 
  // is there to make sure the decay functions work correctly
  hnl.kmom = TLorentzVector(TVector3(0, 0, 0), Constants::Instance().kplus_mass);
  hnl.sec = TLorentzVector(-P * TVector3(1, 0, 0), sqrt(P*P + lep_mass*lep_mass));

  double width = 0.;
  for (const HNLMakeDecay::HNLDecayFunction F: fSelectedDecays) {
    width += ((*this.*F)(hnl)).width;
  }

  double lifetime_ns = Constants::Instance().hbar / width;
  float mean_dist = lifetime_ns * hnl.mom.Gamma() * hnl.mom.Beta() * Constants::Instance().c_cm_per_ns;

  return mean_dist;
}

double HNLMakeDecay::CalculateMaxWeight() {
  double ue4 = fReferenceUE4;
  double um4 = fReferenceUM4;
  double hnl_mass = fReferenceHNLMass;
  double length = fReferenceRayLength;
  double E = fReferenceHNLEnergy;
  double P = sqrt(E*E - hnl_mass * hnl_mass);

  // Compute the boost to get the HNL with the correct energy
  double lep_mass = (um4 > 0) ? Constants::Instance().muon_mass : Constants::Instance().elec_mass;
  double p0 = twobody_momentum(Constants::Instance().kplus_mass, lep_mass, hnl_mass);
  double e0 = sqrt(p0*p0 + hnl_mass*hnl_mass);
  double beta = (-e0*p0 + E*P) / (E*E + p0*p0); 
  TVector3 boost(beta, 0, 0);

  // make a fake HNL flux with these parameters
  MeVPrtlFlux hnl;
  hnl.C1 = ue4;
  hnl.C2 = um4;
  hnl.mass = hnl_mass;
  hnl.mom = TLorentzVector(p0 * TVector3(1, 0, 0), e0);
  hnl.mom.Boost(boost);
  hnl.secondary_pdg = -1; // doesn't affect the end width, just make it viable

  std::cout << "Reference Energy: " << E << std::endl;
  std::cout << "HNL 4P: " << hnl.mom.E() << " " << hnl.mom.Px() << std::endl;

  // Mock up the secondary momenta -- this shouldn't affect width and just 
  // is there to make sure the decay functions work correctly
  hnl.kmom = TLorentzVector(TVector3(0, 0, 0), Constants::Instance().kplus_mass);
  hnl.kmom.Boost(boost);
  hnl.sec = TLorentzVector(-p0 * TVector3(1, 0, 0), sqrt(p0*p0 + lep_mass*lep_mass));
  hnl.sec.Boost(boost);

  double width = 0.;
  for (const HNLMakeDecay::HNLDecayFunction F: fSelectedDecays) {
    width += ((*this.*F)(hnl)).width;
  }

  std::cout << "REFERENCE WIDTH: " << width << std::endl;

  double lifetime_ns = Constants::Instance().hbar / width;
  float mean_dist = lifetime_ns * hnl.mom.Gamma() * hnl.mom.Beta() * Constants::Instance().c_cm_per_ns;

  std::cout << "REFERENCE DECAY LENGTH: " << mean_dist << std::endl;

  return length / mean_dist; 
}


HNLMakeDecay::HNLMakeDecay(fhicl::ParameterSet const &pset):
  IMeVPrtlStage("HNLMakeDecay") 
{
    this->configure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

HNLMakeDecay::~HNLMakeDecay()
{
  gsl_integration_workspace_free(fIntegrator);
}

//------------------------------------------------------------------------------------------------------------------------------------------
void HNLMakeDecay::configure(fhicl::ParameterSet const &pset)
{
  fIntegratorSize = 1000;

  fIntegrator = gsl_integration_workspace_alloc(fIntegratorSize);

  // Setup available decays
  fAvailableDecays["mu_pi"] = &HNLMakeDecay::MuPi;
  fAvailableDecayMasses["mu_pi"] = Constants::Instance().muon_mass + Constants::Instance().piplus_mass;
  fAvailableDecays["e_pi"] = &HNLMakeDecay::EPi;
  fAvailableDecayMasses["e_pi"] = Constants::Instance().elec_mass + Constants::Instance().piplus_mass;
  // TODO: make available
  // fAvailableDecays["nu_mu_mu"] = &HNLMakeDecay::NuMupMum;

  // Select which ones are configued
  fDecayConfig = pset.get<std::vector<std::string>>("Decays");

  for (const std::string &d: fDecayConfig) {
    if (fAvailableDecays.count(d)) {
      fSelectedDecays.push_back(fAvailableDecays.at(d));
    }
    else {
      std::cerr << "ERROR: Selected unavailable decay (" << d << ")" << std::endl;
    }
  }
  fReferenceUE4 = pset.get<double>("ReferenceUE4");
  fReferenceUM4 = pset.get<double>("ReferenceUM4");
  fReferenceHNLMass = pset.get<double>("ReferenceHNLMass");
  fReferenceRayLength = pset.get<double>("ReferenceRayLength");

  fReferenceHNLEnergy = pset.get<float>("ReferenceHNLEnergy", -1);
  fReferenceHNLKaonEnergy = pset.get<float>("ReferenceHNLEnergyFromKaonEnergy", -1.);
  if (fReferenceHNLEnergy < 0. && fReferenceHNLKaonEnergy > 0.) {
    double lep_mass = (fReferenceUE4 > 0) ? Constants::Instance().elec_mass : Constants::Instance().muon_mass;
    fReferenceHNLEnergy = forwardPrtlEnergy(Constants::Instance().kplus_mass, lep_mass, fReferenceHNLMass, fReferenceHNLKaonEnergy);
  }

  fMinDetectorDistance = pset.get<double>("MinDetectorDistance", 100e2); // 100m for NuMI -> SBN/ICARUS

  fMajorana = pset.get<bool>("Majorana");

  fMaxWeight = CalculateMaxWeight();

  // Guard against bad couplings that will mess with decay weighting approximations
  double kdar_decay_length = CalculateKDARDecayLength();
  if (kdar_decay_length < fMinDetectorDistance) {
    throw cet::exception("HNLMakeDecay Tool: BAD CONFIG. Existing configuration results in a KDAR-flux partial decay length (" +
          std::to_string(kdar_decay_length) + ") that is smaller than the configured min distance to the detector (" + std::to_string(fMinDetectorDistance) + ").");
  }
}

bool HNLMakeDecay::Decay(const MeVPrtlFlux &flux, const TVector3 &in, const TVector3 &out, MeVPrtlDecay &decay, double &weight) {
  // Check that the mass/decay configuration is allowed
  bool has_allowed_decay = false;
  for (const std::string &d: fDecayConfig) {
    if (fReferenceHNLMass > fAvailableDecayMasses[d]) {
      has_allowed_decay = true;
      break;
    }
  }

  if (!has_allowed_decay) {
    throw cet::exception("HNLMakeDecay Tool: BAD MASS. Configured mass (" + std::to_string(flux.mass) +
         ") is smaller than any configured decay.");
  }

  // Run the selected decay channels
  std::vector<HNLMakeDecay::DecayFinalState> decays;
  double total_width = 0.;
  for (const HNLMakeDecay::HNLDecayFunction F: fSelectedDecays) {
    decays.push_back((*this.*F)(flux));
    total_width += decays.back().width;
  }

  std::cout << "TOTAL DECAY WIDTH: " << total_width << std::endl;
  if (total_width == 0.) return false;

  // pick one
  double sum_width = 0.;
  int idecay = decays.size()-1;
  double rand = GetRandom();
  for (unsigned i = 0; i < decays.size()-1; i++) {
    sum_width += decays[i].width;
    if (rand < sum_width / total_width) {
      idecay = i;
      break;
    }
  }

  // Pick a random decay position -- assumes it is uniform across the detector
  TVector3 decay_pos = in + GetRandom() * (out - in);

  // Get the decay probability

  // total lifetime
  double lifetime_ns = Constants::Instance().hbar / total_width;

  // multiply by gamma*v to get the length
  float mean_dist = lifetime_ns * flux.mom.Gamma() * flux.mom.Beta() * Constants::Instance().c_cm_per_ns;

  // Get the weight (NOTE: this negelects the probability that the HNL decays before the detector)
  // I.e. it is only valid in the limit mean_dist >> 100m (distance from beam to SBN)
  if (mean_dist < fMinDetectorDistance) {
    std::cerr << "ERROR: bad mean_dist value (" << mean_dist << "). Decay weighting approximations invalid. Ignoring event." << std::endl;
    std::cout << "ERROR: bad mean_dist value (" << mean_dist << "). Decay weighting approximations invalid. Ignoring event." << std::endl;
    return false;
  }

  // saves the weight
  weight = (out - in).Mag() / mean_dist;

  // Save the decay info
  decay.pos = TLorentzVector(decay_pos, TimeOfFlight(flux, decay_pos));
  for (const TLorentzVector &p: decays[idecay].mom) {
    decay.daughter_mom.push_back(p.Vect());
    decay.daughter_e.push_back(p.E());
  }
  decay.daughter_pdg = decays[idecay].pdg;
  decay.decay_width = total_width;
  decay.mean_lifetime = lifetime_ns;
  decay.mean_distance = mean_dist;

  return true;
}

DEFINE_ART_CLASS_TOOL(HNLMakeDecay)

} // namespace ldm
} // namespace evgen
