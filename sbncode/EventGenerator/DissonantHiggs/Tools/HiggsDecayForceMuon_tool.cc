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
#include "IHiggsDecay.h"
#include "../Products/HiggsFlux.h"

// LArSoft includes
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/JamesRandom.h"
#include "CLHEP/Random/RandFlat.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// std includes
#include <string>
#include <iostream>
#include <memory>
#include <utility>

// constants
#define MUON_MASS 105.66 // MeV
#define C_CM_PER_NS 29.9792 

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace evgen {
namespace ldm {
/**
 *  @brief  HiggsDecayForceMuon class definiton
 */
class HiggsDecayForceMuon : public IHiggsDecay {
public:
    /**
     *  @brief  Constructor
     */
    HiggsDecayForceMuon(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~HiggsDecayForceMuon();

    void configure(fhicl::ParameterSet const &pset) override;

    bool Decay(const HiggsFlux &flux, const TVector3 &in, const TVector3 &out, simb::MCTruth &truth, float &weight) override;

    // no weights
    float ConstantWeight() override { return 1.; }
    float MaxWeight() override { return 1.; }

private:
  float fReferenceRayLength;
};

HiggsDecayForceMuon::HiggsDecayForceMuon(fhicl::ParameterSet const &pset):
  IHiggsStage("HiggsDecayForceMuon") 
{
    this->configure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

HiggsDecayForceMuon::~HiggsDecayForceMuon()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
void HiggsDecayForceMuon::configure(fhicl::ParameterSet const &pset)
{
  fReferenceRayLength = pset.get<float>("ReferenceRayLength", -1);
}

// converts a random number (x) between 0 and 1 to a number
// from an exponential distribution with mean forced to lie 
// between a and b
float flat_to_exp_rand(float x, float mean, float a, float b) {
  float A = (1. - exp(-(b-a)/mean));
  return - mean * log(1 - x * A) + a;
}

// returns the weight associated with forcing the decay to happen within a center length
float forcedecay_weight(float mean, float a, float b) {
    return exp(-a/mean) - exp(-b/mean);
}



bool HiggsDecayForceMuon::Decay(const HiggsFlux &flux, const TVector3 &in, const TVector3 &out, simb::MCTruth &truth, float &weight) {

  // get the mean length travelled by the higgs
  float lifetime_ns = 100.; // TODO: set
  // multiply by gamma*v to get the length
  float mean_dist = lifetime_ns * flux.mom.Gamma() * flux.mom.Beta() * C_CM_PER_NS;

  // distance inside detector
  float in_dist = (flux.pos.Vect() - in).Mag();
  float out_dist = (flux.pos.Vect() - out).Mag();

  // compute the decay weight
  // weight = forcedecay_weight(mean_dist, in_dist, out_dist);
  weight = 1.;

  // get the decay location
  float flat_rand = CLHEP::RandFlat::shoot(fEngine, 0, 1.);
  float decay_rand = flat_to_exp_rand(flat_rand, mean_dist, in_dist, out_dist);
  TVector3 decay_pos3 = flux.pos.Vect() + decay_rand * (in - flux.pos.Vect()).Unit();

  // decay time
  float decay_time = flux.pos.T() + decay_rand / (flux.mom.Beta() * C_CM_PER_NS);
  TLorentzVector decay_pos(decay_pos3, decay_time);

  // force the decay to muons
  assert(flux.mom.M() / 2. > MUON_MASS);

  // muon mom+energy in the parent rest-frame
  float muonE_HRF = flux.mom.M() / 2.;
  float muonP_HRF = sqrt(muonE_HRF * muonE_HRF - MUON_MASS * MUON_MASS);

  // make the two muons in the higgs rest-frame with a random angle
  float theta = CLHEP::RandFlat::shoot(fEngine, 0, M_PI);
  float phi = CLHEP::RandFlat::shoot(fEngine, 0, 2*M_PI);

  TVector3 pA_HRF;
  pA_HRF.SetMagThetaPhi(muonP_HRF, theta, phi);

  TVector3 pB_HRF = -pA_HRF;

  TLorentzVector p4A = TLorentzVector(pA_HRF, muonE_HRF);
  TLorentzVector p4B = TLorentzVector(pB_HRF, muonE_HRF);

  // Boost
  p4A.Boost(flux.mom.BoostVector());
  p4B.Boost(flux.mom.BoostVector());

  // make the two MC particles
  simb::MCParticle A(
    0, // track ID
    13, // track pdg
    "primary", // process
    -1, // mother is not available
    MUON_MASS); 

  // add the start point
  A.AddTrajectoryPoint(decay_pos, p4A);

  simb::MCParticle B(
    0, // track ID
    -13, // track pdg
    "primary", // process
    -1, // mother is not available
    MUON_MASS); 

  // add the start point
  B.AddTrajectoryPoint(decay_pos, p4B);

  truth.Add(A);
  truth.Add(B);

  return true;
}

DEFINE_ART_CLASS_TOOL(HiggsDecayForceMuon)

#undef MUON_MASS
#undef C_CM_PER_NS

} // namespace ldm
} // namespace evgen
