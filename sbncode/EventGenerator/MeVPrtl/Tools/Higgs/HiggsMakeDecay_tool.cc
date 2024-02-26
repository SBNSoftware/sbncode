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

// constants
#include "TDatabasePDG.h"

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace evgen {
namespace ldm {
/**
 *  @brief  HiggsMakeDecay class definiton
 *
 *  Implementation of model taken from:
 *      https://arxiv.org/abs/1909.11670
 */
class HiggsMakeDecay : public IMeVPrtlDecay {
public:
    /**
     *  @brief  Constructor
     */
    HiggsMakeDecay(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~HiggsMakeDecay();

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
  double fReferenceHiggsMass;
  double fReferenceHiggsMixing;
  double fReferenceHiggsEnergy;
  double fReferenceHiggsKaonEnergy;

  double fMaxWeight;

  bool fAllowElectronDecay;
  bool fAllowMuonDecay;
  bool fAllowPionDecay;
  bool fAllowPi0Decay;
  bool fAddTimeOfFlight;
  
};

// converts a random number (x) between 0 and 1 to a number
// from an exponential distribution with mean forced to lie 
// between a and b
double flat_to_exp_rand(double x, double mean, double a, double b) {
  double A = (1. - exp(-(b-a)/mean));
  return - mean * log(1 - x * A) + a;
}

// returns the weight associated with forcing the decay to happen within a center length
double forcedecay_weight(double mean, double a, double b) {
    return exp(-a/mean) - exp(-b/mean);
}

double higgs_momentum(double kaon_mass, double pion_mass, double higs_mass) {
  return sqrt(kaon_mass * kaon_mass * kaon_mass * kaon_mass
    -2 * kaon_mass * kaon_mass * pion_mass * pion_mass
    -2 * kaon_mass * kaon_mass * higs_mass * higs_mass
       + pion_mass * pion_mass * pion_mass * pion_mass
       + higs_mass * higs_mass * higs_mass * higs_mass
    -2 * pion_mass * pion_mass * higs_mass * higs_mass) / ( 2 * kaon_mass );
}


// Get the partial width for lepton decays
double LeptonPartialWidth(double lep_mass, double higs_mass, double mixing) {
  if (lep_mass * 2 >= higs_mass) return 0.;

  double higgs_vev = Constants::Instance().higgs_vev;

  double width = (mixing * mixing * lep_mass * lep_mass * higs_mass / (8 * M_PI * higgs_vev * higgs_vev)) * pow(1 - 4 * lep_mass * lep_mass / (higs_mass * higs_mass), 3. / 2.);
  return width;
}

double ElectronPartialWidth(double higs_mass, double mixing) {
  return LeptonPartialWidth(Constants::Instance().elec_mass, higs_mass, mixing);
}

double MuonPartialWidth(double higs_mass, double mixing) {
  return LeptonPartialWidth(Constants::Instance().muon_mass, higs_mass, mixing);
}

double PionPartialWidth(double pion_mass, double higs_mass, double mixing) {
  if (pion_mass * 2 >= higs_mass) return 0.;

  double higgs_vev = Constants::Instance().higgs_vev;

  double form_factor = (2. / 9.) * higs_mass * higs_mass + (11. / 9.) * pion_mass * pion_mass;

  double width = (mixing * mixing * 3 * form_factor * form_factor / (32 * M_PI * higgs_vev * higgs_vev * higs_mass)) * pow(1- 4. * pion_mass * pion_mass / (higs_mass * higs_mass), 1. / 2.);

  return width;
}

double PiPlusPartialWidth(double higs_mass, double mixing) {
  return PionPartialWidth(Constants::Instance().piplus_mass, higs_mass, mixing);
}

double PiZeroPartialWidth(double higs_mass, double mixing) {
  return PionPartialWidth(Constants::Instance().pizero_mass, higs_mass, mixing);
}

HiggsMakeDecay::HiggsMakeDecay(fhicl::ParameterSet const &pset):
  IMeVPrtlStage("HiggsMakeDecay") 
{
    this->configure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

HiggsMakeDecay::~HiggsMakeDecay()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
void HiggsMakeDecay::configure(fhicl::ParameterSet const &pset)
{
  fReferenceRayLength = pset.get<double>("ReferenceRayLength", -1);
  fReferenceRayDistance = pset.get<double>("ReferenceRayDistance", 0.);
  fReferenceHiggsMass = pset.get<double>("ReferenceHiggsMass", -1);
  fReferenceHiggsMixing = pset.get<double>("ReferenceHiggsMixing", -1);
  fReferenceHiggsEnergy = pset.get<double>("ReferenceHiggsEnergy", -1);
  fReferenceHiggsKaonEnergy = pset.get<double>("ReferenceHiggsEnergyFromKaonEnergy", -1.);

  fAllowElectronDecay = pset.get<bool>("AllowElectronDecay", true);
  fAllowMuonDecay = pset.get<bool>("AllowMuonDecay", true);
  fAllowPionDecay = pset.get<bool>("AllowPionDecay", true);
  fAllowPi0Decay = pset.get<bool>("AllowPi0Decay", true);
  fAddTimeOfFlight = pset.get<bool>("AddTimeOfFlight", true);

  if (fReferenceHiggsEnergy < 0. && fReferenceHiggsKaonEnergy > 0.) {
    fReferenceHiggsEnergy = std::min(forwardPrtlEnergy(Constants::Instance().kplus_mass, Constants::Instance().piplus_mass, fReferenceHiggsMass, fReferenceHiggsKaonEnergy),
                                     forwardPrtlEnergy(Constants::Instance().klong_mass, Constants::Instance().pizero_mass, fReferenceHiggsMass, fReferenceHiggsKaonEnergy));
  }

  // if configured to, divide out some of the decay weight
  if (fReferenceRayLength > 0. && fReferenceHiggsMass > 0. && fReferenceHiggsMixing > 0. && fReferenceHiggsEnergy > 0.) {
    // Get each partial width
    double width_elec = ElectronPartialWidth(fReferenceHiggsMass, fReferenceHiggsMixing);
    double width_muon = MuonPartialWidth(fReferenceHiggsMass, fReferenceHiggsMixing);
    double width_piplus = PiPlusPartialWidth(fReferenceHiggsMass, fReferenceHiggsMixing);
    double width_pizero = PiZeroPartialWidth(fReferenceHiggsMass, fReferenceHiggsMixing);

    // total lifetime
    double lifetime_ns = Constants::Instance().hbar / (width_elec + width_muon + width_piplus + width_pizero);

    // multiply by gamma*v to get the length
    double gamma_v = sqrt(fReferenceHiggsEnergy * fReferenceHiggsEnergy - fReferenceHiggsMass * fReferenceHiggsMass) * Constants::Instance().c_cm_per_ns / fReferenceHiggsMass;
    double mean_dist = lifetime_ns * gamma_v;

    // compute the decay weight
    fMaxWeight = forcedecay_weight(mean_dist, fReferenceRayDistance, fReferenceRayDistance + fReferenceRayLength);

    // Scale by the allowed BR
    fMaxWeight *= (width_elec*fAllowElectronDecay + width_muon*fAllowMuonDecay + width_piplus*fAllowPionDecay + width_pizero*fAllowPi0Decay) / (width_elec + width_muon + width_piplus + width_pizero);

  }
  else {
    fMaxWeight = -1.;
  }
}

int HiggsMakeDecay::RandDaughter(double elec_width, double muon_width, double piplus_width, double pizero_width) {
  double total_width = elec_width + muon_width + piplus_width + pizero_width;

  double flat_rand = CLHEP::RandFlat::shoot(fEngine, 0, 1.);

  if (flat_rand < elec_width / total_width) {
    return 11;
  }
  else if (flat_rand < (elec_width + muon_width) / total_width) {
    return 13;
  }
  else if (flat_rand < (elec_width + muon_width + piplus_width) / total_width) {
    return 211;
  }
  else {
    return 111;
  }

  // unreachable
  return -1;
}

bool HiggsMakeDecay::Decay(const MeVPrtlFlux &flux, const TVector3 &in, const TVector3 &out, MeVPrtlDecay &decay, double &weight) {
  // Handle bad mass value
  if (flux.mass < 2*Constants::Instance().elec_mass) {
    throw cet::exception("HiggsMakeDecay Tool: BAD MASS. Configured mass (" + std::to_string(flux.mass) +
         ") is smaller than lowest mass available decay e+e- (" + std::to_string(2*Constants::Instance().elec_mass) +")");
  } 

  double mixing = flux.C1; // constant-1 saves the mixing

  // Get each partial width
  double width_elec = ElectronPartialWidth(flux.mass, mixing);
  double width_muon = MuonPartialWidth(flux.mass, mixing);
  double width_piplus = PiPlusPartialWidth(flux.mass, mixing);
  double width_pizero = PiZeroPartialWidth(flux.mass, mixing);

  // total lifetime
  double lifetime_ns = Constants::Instance().hbar / (width_elec + width_muon + width_piplus + width_pizero);

  // multiply by gamma*v to get the length
  double mean_dist = lifetime_ns * flux.mom.Gamma() * flux.mom.Beta() * Constants::Instance().c_cm_per_ns;

  // distance inside detector
  double in_dist = (flux.pos.Vect() - in).Mag();
  double out_dist = (flux.pos.Vect() - out).Mag();

  // compute the decay weight
  weight = forcedecay_weight(mean_dist, in_dist, out_dist);

  // if the weight is literally zero, then there is no way in hell
  // that this higgs is going to decay in the detector -- reject it
  if (weight == 0) {
    return false;
  }

  double partial_to_total = (width_elec*fAllowElectronDecay + width_muon*fAllowMuonDecay + width_piplus*fAllowPionDecay + width_pizero*fAllowPi0Decay) / (width_elec + width_muon + width_piplus + width_pizero);

  // Scale by the allowed BR
  weight *= partial_to_total;

  // get the decay location
  double flat_rand = CLHEP::RandFlat::shoot(fEngine, 0, 1.);
  double decay_rand = flat_to_exp_rand(flat_rand, mean_dist, in_dist, out_dist);
  TVector3 decay_pos3 = flux.pos.Vect() + decay_rand * (in - flux.pos.Vect()).Unit();

  // decay time
  double decay_time = fAddTimeOfFlight ? TimeOfFlight(flux, decay_pos3) : flux.pos.T() /* no TOF: set time equal to parent */;
  TLorentzVector decay_pos(decay_pos3, decay_time);

  // get the decay type
  int daughter_pdg = RandDaughter(width_elec*fAllowElectronDecay, width_muon*fAllowMuonDecay, width_piplus*fAllowPionDecay, width_pizero*fAllowPi0Decay);

  double daughter_mass = TDatabasePDG::Instance()->GetParticle(daughter_pdg)->Mass();

  // daughter mom+energy in the parent rest-frame
  double daughterE_HRF = flux.mass / 2.;
  double daughterP_HRF = sqrt(daughterE_HRF * daughterE_HRF - daughter_mass * daughter_mass);

  // make the two daughters in the higgs rest-frame with a random direction
  TVector3 pA_HRF = RandomUnitVector() * daughterP_HRF;

  TVector3 pB_HRF = -pA_HRF;

  TLorentzVector p4A = TLorentzVector(pA_HRF, daughterE_HRF);
  TLorentzVector p4B = TLorentzVector(pB_HRF, daughterE_HRF);

  // Boost
  p4A.Boost(flux.mom.BoostVector());
  p4B.Boost(flux.mom.BoostVector());

  // save the decay info
  decay.total_decay_width = width_elec + width_muon + width_piplus + width_pizero; 
  decay.total_mean_lifetime = lifetime_ns;
  decay.total_mean_distance = mean_dist; 
  decay.allowed_decay_fraction = partial_to_total;

  decay.pos = decay_pos;

  decay.daughter_mom.push_back(p4A.Vect());
  decay.daughter_e.push_back(p4A.E());
  decay.daughter_pdg.push_back(daughter_pdg);

  decay.daughter_mom.push_back(p4B.Vect());
  decay.daughter_e.push_back(p4B.E());
  // daughter B is anti-particle 
  if (daughter_pdg == 111) { // pi0 is its own anti-particle
    decay.daughter_pdg.push_back(daughter_pdg);
  }
  else {
    decay.daughter_pdg.push_back(-daughter_pdg);
  }

  return true;
}

DEFINE_ART_CLASS_TOOL(HiggsMakeDecay)

} // namespace ldm
} // namespace evgen
