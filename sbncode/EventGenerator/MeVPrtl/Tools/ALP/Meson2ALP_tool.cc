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
#include "CLHEP/Random/RandFlat.h"

// local includes
#include "sbncode/EventGenerator/MeVPrtl/Tools/Constants.h"
#include "sbncode/EventGenerator/MeVPrtl/Tools/IMeVPrtlFlux.h"

#include "sbnobj/Common/EventGen/MeVPrtl/MeVPrtlTruth.h"
#include "sbnobj/Common/EventGen/MeVPrtl/MeVPrtlFlux.h"
#include "sbnobj/Common/EventGen/MeVPrtl/MesonParent.h"

// LArSoft includes
#include "nugen/EventGeneratorBase/GENIE/EvtTimeShiftFactory.h"
#include "nugen/EventGeneratorBase/GENIE/EvtTimeShiftI.h"

// ROOT
#include "TVector3.h"

// std includes
#include <string>
#include <iostream>
#include <memory>

#include "TDatabasePDG.h"

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace evgen {
namespace ldm {
/**
 *  @brief  Meson2ALP class definiton
 *  Implementation of model taken from:
 *      https://arxiv.org/abs/1909.11670
 */
class Meson2ALP : public IMeVPrtlFlux
{
public:
    /**
     *  @brief  Constructor
     */
    Meson2ALP(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~Meson2ALP();

    bool MakeFlux(const simb::MCFlux &flux, MeVPrtlFlux &alp, double &weight) override;
    void configure(const fhicl::ParameterSet&) override;

    double MaxWeight() override; 

    double EtaBR() const;
    double EtaPBR() const;
    double Pi0BR() const;

private:
  // config
  double fM; //!< Mass of axion [GeV]
  double ffa; //!< Axion decay constant [GeV]
  double fcAl; //!< Axial coupling to leptons
  double fcG; //!< Axion coupling to gluon
  double fcB; //!< Axion coupling to U(1) B boson (before EW symmetry breaking)
  double fcW; //!< Axion coupling to SU(2) W bosons (before EW symmetry breaking)

  // branching ratios
  double fPi0BR;
  double fEtaBR;
  double fEtaPBR;

  double fMaxWeightPi0;
  double fMaxWeightEta;
  double fMaxWeightEtaP;

  // weight
  double fMaxWeight;
};

Meson2ALP::Meson2ALP(fhicl::ParameterSet const &pset):
  IMeVPrtlStage("Meson2ALP"), 
  IMeVPrtlFlux(pset) 
{
    this->configure(pset);

}

//------------------------------------------------------------------------------------------------------------------------------------------

Meson2ALP::~Meson2ALP() {}

//------------------------------------------------------------------------------------------------------------------------------------------
void Meson2ALP::configure(fhicl::ParameterSet const &pset)
{
  fM = pset.get<double>("M");
  ffa = pset.get<double>("fa");
  fcAl = pset.get<double>("cAl");
  fcG = pset.get<double>("cG");
  fcB = pset.get<double>("cB");
  fcW = pset.get<double>("cW");

  fMaxWeightPi0 = pset.get<double>("MaxWeightPi0", 1.);
  fMaxWeightEta = pset.get<double>("MaxWeightEta", 1.);
  fMaxWeightEtaP = pset.get<double>("MaxWeightEtaP", 1.);

  fPi0BR = Pi0BR();
  fEtaBR = EtaBR();
  fEtaPBR = EtaPBR();

  fMaxWeight = std::max(fPi0BR*fMaxWeightPi0, std::max(fEtaBR*fMaxWeightEta, fEtaPBR*fMaxWeightEtaP));

  std::cout << "Branching Ratios -- " << "Pi0: " << fPi0BR << " Eta: " << fEtaBR << " Eta': " << fEtaPBR << std::endl;
}

double Meson2ALP::MaxWeight() {
  return fMaxWeight;
}

double Meson2ALP::EtaBR() const {
  double eta_mass = Constants::Instance().eta_mass;
  double pzero_mass = Constants::Instance().pizero_mass;
  double fpion = Constants::Instance().fpion / sqrt(2); // fpi ~ 93MeV convention

  double mixing_angle = (1. / sqrt(6)) * (fcG * fpion / ffa) * (fM*fM - (4./9.)*pzero_mass*pzero_mass) / (fM*fM - eta_mass*eta_mass);
  double qcd_rate_f = 1.;
  if (fM > eta_mass) qcd_rate_f = pow(fM/eta_mass, -1.6);

  return mixing_angle*mixing_angle*qcd_rate_f;
}

double Meson2ALP::EtaPBR() const {
  double etap_mass = Constants::Instance().etap_mass;
  double pzero_mass = Constants::Instance().pizero_mass;
  double fpion = Constants::Instance().fpion / sqrt(2); // fpi ~ 93MeV convention

  double mixing_angle = (1. / sqrt(12)) * (fcG * fpion / ffa) * (fM*fM - (16./9.)*pzero_mass*pzero_mass) / (fM*fM - etap_mass*etap_mass);
  double qcd_rate_f = 1.;
  if (fM > etap_mass) qcd_rate_f = pow(fM/etap_mass, -1.6);

  return mixing_angle*mixing_angle*qcd_rate_f;
}

double Meson2ALP::Pi0BR() const {
  double pzero_mass = Constants::Instance().pizero_mass;
  double fpion = Constants::Instance().fpion / sqrt(2); // fpi ~ 93MeV convention

  double mixing_angle = (1./6) * (fcG * fpion / ffa) * fM*fM / (fM*fM - pzero_mass*pzero_mass);

  double qcd_rate_f = 1.;
  if (fM > pzero_mass) qcd_rate_f = pow(fM/pzero_mass, -1.6);

  return mixing_angle*mixing_angle*qcd_rate_f;
}

bool Meson2ALP::MakeFlux(const simb::MCFlux &flux, evgen::ldm::MeVPrtlFlux &alp, double &weight) {
  evgen::ldm::MesonParent meson(flux);

  TLorentzVector Beam4 = BeamOrigin();
  // get position in detector frame
  alp.pos_beamcoord = meson.pos;
  alp.pos = meson.pos;
  alp.pos.Transform(fBeam2Det);
  alp.pos += Beam4;

  // Energy is same as for meson (don't worry about momentum conservation)
  double alp_energy = meson.mom.E();
  double alp_momentum = sqrt(alp_energy*alp_energy - fM*fM);

  // Momentum 4-vector
  TLorentzVector mom(meson.mom.Vect().Unit()*alp_momentum, alp_energy);

  alp.mom_beamcoord = mom;
  // rotate to detector frame
  alp.mom = mom;
  alp.mom.Transform(fBeam2Det);

  // These are interactions, not decays, 
  // so the "parent" here is a proton in the z-direction
  alp.mmom_beamcoord = TLorentzVector(TVector3(0, 0, 120), sqrt(120*120 + 1));
  alp.mmom = alp.mmom_beamcoord;
  alp.mmom.Transform(fBeam2Det);

  // The weight is the importance weight times the branching-ratio weight 
  weight = meson.weight;
  if (meson.meson_pdg == 111 /* pi0 */) {
    weight = weight * fPi0BR;
  }
  else if (meson.meson_pdg == 221 /* eta */) {
    weight = weight * fEtaBR;
  }
  else if (meson.meson_pdg == 331 /* eta' */) { 
    weight = weight * fEtaPBR;
  }
  else return false;

  // set the mixing
  alp.C1 = ffa;
  alp.C2 = fcAl;
  alp.C3 = fcG;
  alp.C4 = fcB;
  alp.C5 = fcW;
  alp.mass = fM;

  alp.meson_pdg = meson.meson_pdg;
  alp.generator = evgen::ldm::kALP;

  // No secondary
  alp.secondary_pdg = -1;
  // No equivalent neutrino
  alp.equiv_enu = -1;

  return true;
}

DEFINE_ART_CLASS_TOOL(Meson2ALP)

} // namespace ldm
} // namespace evgen
