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
#include "sbncode/EventGenerator/MeVPrtl/Tools/IMeVPrtlFlux.h"
#include "sbncode/EventGenerator/MeVPrtl/Tools/Constants.h"

#include "sbnobj/Common/EventGen/MeVPrtl/MeVPrtlFlux.h"
#include "sbnobj/Common/EventGen/MeVPrtl/MesonParent.h"

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
 *  @brief  Tau2HNLFlux class definiton
 *
 *  Implementation of Kaon->HNL branching ratio taken from:
 *      arXiv:1912.07622
 */
class Tau2HNLFlux : public IMeVPrtlFlux
{
public:
    /**
     *  @brief  Constructor
     */
    Tau2HNLFlux(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~Tau2HNLFlux();

    bool MakeFlux(const simb::MCFlux &flux, MeVPrtlFlux &hnl, double &weight) override;
    void configure(const fhicl::ParameterSet&) override;

    double MaxWeight() override;

private:
  // config
  double fM; //!< Mass of HNL [GeV]
  double fMagUt4;
};

Tau2HNLFlux::Tau2HNLFlux(fhicl::ParameterSet const &pset):
  IMeVPrtlStage("Tau2HNLFlux"), 
  IMeVPrtlFlux(pset) 
{
    this->configure(pset);

}

//------------------------------------------------------------------------------------------------------------------------------------------

Tau2HNLFlux::~Tau2HNLFlux()
{
}

// helper functions
//
// tau -> pi + HNL
double hnl_momentum(double tau_mass, double pion_mass, double hnl_mass) {
  if (tau_mass - pion_mass < hnl_mass) return -1.;

  return sqrt(tau_mass * tau_mass * tau_mass * tau_mass 
    -2 * tau_mass * tau_mass * pion_mass * pion_mass
    -2 * tau_mass * tau_mass * hnl_mass * hnl_mass
       + pion_mass * pion_mass * pion_mass * pion_mass 
       + hnl_mass * hnl_mass * hnl_mass * hnl_mass 
    -2 * pion_mass * pion_mass * hnl_mass * hnl_mass) / ( 2 * tau_mass );
}

/* Polarization of a HNL emerging from a two body meson decay (M-> l N)
   m_l is the mass of the lepton and m_M the mass of the parent meson
   See arXiv:2109.10358 for more details.
   Due to crossing symmetry the polarization of the tau -> m N gives the same result */
double PolHNL(double m_HNL, double m_l, double m_M)
{
  double yl = m_l / m_M;
  double yHNL =  m_HNL / m_M;
  double numterm = (yl*yl - yHNL * yHNL) * sqrt(pow(yHNL,4) + (1.0 - yl * yl)*(1.0 - yl * yl) - 2.0 * yHNL * yHNL * (1.0 + yl * yl));
  double denterm = pow(yl,4) + pow(yHNL,4) - 2.0 * yl * yl * yHNL * yHNL - yl * yl - yHNL * yHNL;
  return numterm / denterm;
}

double BranchingRatio(double hnl_mass, double u4) {
  double Gfermi = Constants::Instance().Gfermi;
  double Vud2 = Constants::Instance().abs_Vud_squared;
  double fpion = Constants::Instance().fpion;
  double tau_mass = Constants::Instance().tau_mass;
  double tau_lifetime = Constants::Instance().tau_lifetime;
  double pi_mass = Constants::Instance().piplus_mass;
  double hbar = Constants::Instance().hbar;
  
  double coupling = 0.9*tau_lifetime*u4*Gfermi*Gfermi*Vud2*fpion*fpion*tau_mass*tau_mass*tau_mass / (16*M_PI) / hbar; 

  double hnl_mass2 = hnl_mass*hnl_mass;
  double tau_mass2 = tau_mass*tau_mass;
  double pi_mass2 = pi_mass*pi_mass;

  double kinematic_factor = \
    ((1 - hnl_mass2 / tau_mass2) * (1 - hnl_mass2 / tau_mass2) \
     - (pi_mass2 / tau_mass2) * (1 + hnl_mass2 / tau_mass2)) \
    *sqrt((1 - (pi_mass - hnl_mass)*(pi_mass - hnl_mass)/tau_mass2)*(1 - (hnl_mass+pi_mass)*(hnl_mass+pi_mass)/tau_mass2));
  
  return coupling * kinematic_factor;
}

//------------------------------------------------------------------------------------------------------------------------------------------
void Tau2HNLFlux::configure(fhicl::ParameterSet const &pset)
{
  fM = pset.get<double>("M");
  fMagUt4 = pset.get<double>("MagUt4");

  double max_mass = Constants::Instance().tau_mass - Constants::Instance().piplus_mass;

  if (fM > max_mass) {
    throw cet::exception("Tau2HNLFlux Tool: BAD MASS. Configured mass (" + std::to_string(fM) +
         ") is larger than maximum allowed by enabled couplings (" + std::to_string(max_mass) +  ").");
  }

}

double Tau2HNLFlux::MaxWeight() { 
  return BranchingRatio(fM, fMagUt4);
}

bool Tau2HNLFlux::MakeFlux(const simb::MCFlux &flux, evgen::ldm::MeVPrtlFlux &hnl, double &weight) {
  // make the tau parent
  evgen::ldm::MesonParent tau(flux);
  if (abs(tau.meson_pdg) != 15) return false; // Only take taus

  TLorentzVector Beam4 = BeamOrigin();

  // get position in detector frame
  hnl.pos_beamcoord = tau.pos;
  hnl.pos = tau.pos;
  hnl.pos.Transform(fBeam2Det);
  hnl.pos += Beam4;

  double br = BranchingRatio(fM, fMagUt4);

  std::cout << "BR: " << br << std::endl;

  // ignore if we can't make this hnl
  // Ignore if branching ratio is exactly 0.
  if (br == 0.) return false;

  // get the momentum direction in the tau parent rest frame
  double p = hnl_momentum(Constants::Instance().tau_mass, Constants::Instance().piplus_mass, fM);
  double e = sqrt(p*p + fM * fM);
  // Two-body decays are isotropic
  hnl.mom = TLorentzVector(p*RandomUnitVector(), e);

  // boost to lab frame
  TLorentzVector mom = hnl.mom;
  mom.Boost(tau.mom.BoostVector());

  hnl.mom_beamcoord = mom;
  // rotate to detector frame
  hnl.mom = mom;
  hnl.mom.Transform(fBeam2Det);

  hnl.mmom_beamcoord = tau.mom;
  // also save the tau momentum in the detector frame
  hnl.mmom = tau.mom;
  hnl.mmom.Transform(fBeam2Det);

  // and save the secondary momentum
  hnl.sec = hnl.mmom - hnl.mom;
  hnl.sec_beamcoord = hnl.mmom_beamcoord - hnl.mom_beamcoord;

  // The weight is the importance weight times the branching-ratio weight 
  weight = tau.weight * br;

  // set the mixing
  hnl.C1 = 0.;
  hnl.C2 = 0.;
  hnl.C3 = fMagUt4;
  hnl.mass = fM;

  hnl.meson_pdg = tau.meson_pdg;
  hnl.secondary_pdg = 211 * (tau.meson_pdg > 0 ? 1 : -1);
  hnl.generator = 1; // kHNL

  // equivalent neutrino energy -- doesn't apply here
  hnl.equiv_enu = -1;


// Get the HNL Polarization 
  double meson_mass = Constants::Instance().kplus_mass;
  if(abs(hnl.meson_pdg) == 321) {
    meson_mass = Constants::Instance().kplus_mass;
  }else if(abs(hnl.meson_pdg) == 211) {
    meson_mass = Constants::Instance().piplus_mass;
  }
  

  hnl.polarization = PolHNL(fM , Constants::Instance().tau_mass, meson_mass);

  return true;
}

DEFINE_ART_CLASS_TOOL(Tau2HNLFlux)

} // namespace ldm
} // namespace evgen
