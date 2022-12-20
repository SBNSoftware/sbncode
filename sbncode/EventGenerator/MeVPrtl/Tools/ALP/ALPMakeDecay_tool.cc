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

// local includes
#include "sbncode/EventGenerator/MeVPrtl/Tools/IMeVPrtlDecay.h"
#include "sbncode/EventGenerator/MeVPrtl/Tools/Constants.h"

#include "sbnobj/Common/EventGen/MeVPrtl/MeVPrtlFlux.h"

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
#include <complex>
#include <cmath>

// constants
#include "TDatabasePDG.h"

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace evgen {
namespace ldm {
/**
 *  @brief  ALPMakeDecay class definiton
 *
 *  Implementation of model taken from:
 *      https://arxiv.org/abs/1909.11670
 */
class ALPMakeDecay : public IMeVPrtlDecay {
public:
    /**
     *  @brief  Constructor
     */
    ALPMakeDecay(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~ALPMakeDecay();

    void configure(fhicl::ParameterSet const &pset) override;

    bool Decay(const MeVPrtlFlux &flux, const TVector3 &in, const TVector3 &out, MeVPrtlDecay &decay, double &weight) override;
    int RandDaughter(double elec_width, double muon_width, double piplus_width, double pizero_width);

    // returns the max weight of configured
    double MaxWeight() override { 
      return fMaxWeight; 
    }

private:
  double fReferenceRayLength;
  double fReferenceRayDistance;
  double fReferenceALPMass;
  double fReferenceALPcAl;
  double fReferenceALPDecayConstant;
  double fReferenceALPEnergy;
  bool fAllowEMDecay;

  double fMaxWeight;
};

// converts a random number (x) between 0 and 1 to a number
// from an exponential distribution with mean forced to lie 
// between a and b
double flat_to_exp_rand(double x, double mean, double a, double b) {
  double A = (1. - exp(-(b-a)/mean));
  return - mean * log(1. - x * A) + a;
}

// returns the weight associated with forcing the decay to happen within a center length
double forcedecay_weight(double mean, double a, double b) {
    return exp(-a/mean) - exp(-b/mean);
}

// Get the partial width for lepton decays
// a -> u+u
double MuonPartialWidth(double alp_mass, double cAl, double fA) {
  double muon_mass = Constants::Instance().muon_mass;
  if (muon_mass * 2 >= alp_mass) return 0.;

  double width = (cAl*cAl * alp_mass * muon_mass*muon_mass) / (8 * M_PI * fA*fA)*sqrt(1. - 4. * muon_mass*muon_mass / (alp_mass*alp_mass));
  return width;
}

std::complex<double> B1(double x) {
  std::complex<double> g = 0.;
  if (x >= 1.) {
    g = asin(1./sqrt(x));
  }
  else {
    g = std::complex<double>(M_PI/2., std::log((1+sqrt(1-x)) / (1-sqrt(1-x))));
  }
  return 1. - x*g*g;
}

// a -> gamma+gamma
double GammaPartialWidth(double alp_mass, double c1, double c2, double c3, double cAl, double fa) {
  double pion_mass = Constants::Instance().pizero_mass;
  double eta_mass = Constants::Instance().eta_mass;
  double etap_mass = Constants::Instance().etap_mass;
  double cgluon = c3*(-1.92 /*low-E QCD*/ + 
                     (1./3.)*alp_mass*alp_mass / (alp_mass*alp_mass - pion_mass*pion_mass) /* pion */ +
                     (8./9.)*(alp_mass*alp_mass - (4./9.)*pion_mass*pion_mass) / (alp_mass*alp_mass - eta_mass * eta_mass) /* eta */ +
                     (7./9.)*(alp_mass*alp_mass - (16./9.)*pion_mass*pion_mass) / (alp_mass*alp_mass - etap_mass*etap_mass) /* eta' */); 

  double elec_mass = Constants::Instance().elec_mass;
  double muon_mass = Constants::Instance().muon_mass;
  double tau_mass = Constants::Instance().tau_mass;
  std::complex<double> clep = 2*cAl*(B1(4*elec_mass*elec_mass / (alp_mass*alp_mass)) +
                       B1(4*muon_mass*muon_mass / (alp_mass*alp_mass)) +
                       B1(4*tau_mass*tau_mass / (alp_mass*alp_mass)));

  std::complex<double> cgamma = (5./3.)*c1 + c2 + cgluon + clep;

  double fsc = Constants::Instance().fine_structure_constant;
  double width = fsc*fsc*std::norm(cgamma)*alp_mass*alp_mass*alp_mass/(256*M_PI*M_PI*M_PI*fa*fa);

  return width;
}

// a -> pi+pi-pi0 

// a -> pi0pi0pi0

// a -> pi0pi0pi0

// a-> pi0pi0gamma

// a-> pi+pi-gamma

ALPMakeDecay::ALPMakeDecay(fhicl::ParameterSet const &pset):
  IMeVPrtlStage("ALPMakeDecay") 
{
    this->configure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

ALPMakeDecay::~ALPMakeDecay()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
void ALPMakeDecay::configure(fhicl::ParameterSet const &pset)
{
  fReferenceRayLength = pset.get<double>("ReferenceRayLength", -1);
  fReferenceRayDistance = pset.get<double>("ReferenceRayDistance", 0.);
  fReferenceALPMass = pset.get<double>("ReferenceALPMass", -1);
  fReferenceALPcAl = pset.get<double>("ReferenceALPcAl", -1);
  fReferenceALPDecayConstant = pset.get<double>("ReferenceALPDecayConstant", -1);
  fReferenceALPEnergy = pset.get<double>("ReferenceALPEnergy", -1.);
  fAllowEMDecay = pset.get<bool>("AllowEMDecay");

  // if configured to, divide out some of the decay weight
  if (fReferenceRayLength > 0. && fReferenceALPMass > 0. && fReferenceALPDecayConstant >= 0. && fReferenceALPcAl >= 0. && fReferenceALPEnergy > 0.) {
    // Get each partial width
    double width_muon = MuonPartialWidth(fReferenceALPMass, fReferenceALPcAl, fReferenceALPDecayConstant);
    double width_gamma = GammaPartialWidth(fReferenceALPMass, 1, 1, 1, fReferenceALPcAl, fReferenceALPDecayConstant);

    std::cout << "MU WIDTH: " << width_muon << " GAMMA WIDTH: " << width_gamma << std::endl;
    double partial_width = width_muon + width_gamma*fAllowEMDecay;
    double total_width = width_muon + width_gamma;

    // total lifetime
    double lifetime_ns = Constants::Instance().hbar / total_width;

    // multiply by gamma*v to get the length
    double gamma_v = sqrt(fReferenceALPEnergy * fReferenceALPEnergy - fReferenceALPMass * fReferenceALPMass) * Constants::Instance().c_cm_per_ns / fReferenceALPMass;
    double mean_dist = lifetime_ns * gamma_v;

    // compute the decay weight
    fMaxWeight = forcedecay_weight(mean_dist, fReferenceRayDistance, fReferenceRayDistance + fReferenceRayLength) * partial_width / total_width;
  }
  else {
    std::cout << "ALPMakeDecay: No max weight!! Cannot de-weight.\n";
    fMaxWeight = -1.;
  }
}

bool ALPMakeDecay::Decay(const MeVPrtlFlux &flux, const TVector3 &in, const TVector3 &out, MeVPrtlDecay &decay, double &weight) {
  // Handle bad mass value
  if (flux.mass < 2*Constants::Instance().elec_mass) {
    throw cet::exception("ALPMakeDecay Tool: BAD MASS. Configured mass (" + std::to_string(flux.mass) +
         ") is smaller than lowest mass available decay e+e- (" + std::to_string(2*Constants::Instance().elec_mass) +")");
  } 

  double fA = flux.C1; 
  double cAl = flux.C2;

  // Get each partial width
  double width_muon = MuonPartialWidth(flux.mass, cAl, fA);
  double width_gamma = GammaPartialWidth(flux.mass, 1, 1, 1, cAl, fA);

  double total_width = width_muon + width_gamma;
  double partial_width = width_muon + width_gamma*fAllowEMDecay;

  double partial_to_total = partial_width / total_width;

  // total lifetime
  double lifetime_ns = Constants::Instance().hbar / total_width;

  // multiply by gamma*v to get the length
  double mean_dist = lifetime_ns * flux.mom.Gamma() * flux.mom.Beta() * Constants::Instance().c_cm_per_ns;

  // distance inside detector
  double in_dist = (flux.pos.Vect() - in).Mag();
  double out_dist = (flux.pos.Vect() - out).Mag();

  // compute the decay weight
  // Also weight by the probability of decay to muons
  weight = forcedecay_weight(mean_dist, in_dist, out_dist)*partial_width / total_width;

  // if the weight is literally zero, then there is no way in hell
  // that this alp is going to decay in the detector -- reject it
  if (weight == 0) {
    return false;
  }

  // Scale by the allowed BR

  // get the decay location
  double flat_rand = CLHEP::RandFlat::shoot(fEngine, 0, 1.);
  double decay_rand = flat_to_exp_rand(flat_rand, mean_dist, in_dist, out_dist);
  TVector3 decay_pos3 = flux.pos.Vect() + decay_rand * (in - flux.pos.Vect()).Unit();

  // decay time
  double decay_time = TimeOfFlight(flux, decay_pos3);
  TLorentzVector decay_pos(decay_pos3, decay_time);

  // get the decay type
  bool daughter_is_muon = GetRandom() < (width_muon / partial_width);
  int daughter_pdgA = daughter_is_muon ? 13 : 22;
  int daughter_pdgB = daughter_is_muon ? -13 : 22; // anti-particle of A
  double daughter_mass = daughter_is_muon ? Constants::Instance().muon_mass : 0;

  // daughter mom+energy in the parent rest-frame
  double daughterE_HRF = flux.mass / 2.;
  double daughterP_HRF = sqrt(daughterE_HRF * daughterE_HRF - daughter_mass * daughter_mass);

  // make the two daughters in the alp rest-frame with a random direction
  TVector3 pA_HRF = RandomUnitVector() * daughterP_HRF;

  TVector3 pB_HRF = -pA_HRF;

  TLorentzVector p4A = TLorentzVector(pA_HRF, daughterE_HRF);
  TLorentzVector p4B = TLorentzVector(pB_HRF, daughterE_HRF);

  // Boost
  p4A.Boost(flux.mom.BoostVector());
  p4B.Boost(flux.mom.BoostVector());

  // save the decay info
  decay.total_decay_width = total_width; 
  decay.total_mean_lifetime = lifetime_ns;
  decay.total_mean_distance = mean_dist; 
  decay.allowed_decay_fraction = partial_to_total;

  decay.pos = decay_pos;

  decay.daughter_mom.push_back(p4A.Vect());
  decay.daughter_e.push_back(p4A.E());
  decay.daughter_pdg.push_back(daughter_pdgA);

  decay.daughter_mom.push_back(p4B.Vect());
  decay.daughter_e.push_back(p4B.E());
  decay.daughter_pdg.push_back(daughter_pdgB);

  return true;
}

DEFINE_ART_CLASS_TOOL(ALPMakeDecay)

} // namespace ldm
} // namespace evgen
