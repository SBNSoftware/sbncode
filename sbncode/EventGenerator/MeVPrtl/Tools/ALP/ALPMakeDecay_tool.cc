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
#include "sbncode/EventGenerator/MeVPrtl/Tools/ALP/ThreeBodyIntegrator.h"

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
#include <functional>

#include <gsl/gsl_integration.h>

// constants
#include "TDatabasePDG.h"

using namespace std::complex_literals;

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

    // returns the max weight of configured
    double MaxWeight() override { 
      return fMaxWeight; 
    }

private:
  double fReferenceRayLength;
  double fReferenceRayDistance;
  double fReferenceALPMass;
  double fReferenceALPcAl;
  double fReferenceALPcG;
  double fReferenceALPcB;
  double fReferenceALPcW;
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
    g = std::complex<double>(M_PI/2., std::log((1+sqrt(1-x)) / (1-sqrt(1-x)))/2);
  }
  return 1. - x*g*g;
}

// a -> gamma+gamma
double GammaPartialWidth(double alp_mass, double cB, double cW, double cG, double cAl, double fa) {
  double pion_mass = Constants::Instance().pizero_mass;
  double eta_mass = Constants::Instance().eta_mass;
  double etap_mass = Constants::Instance().etap_mass;
  double cgluon = cG*(-1.92 /*low-E QCD*/ + 
                     (1./3.)*alp_mass*alp_mass / (alp_mass*alp_mass - pion_mass*pion_mass) /* pion */ +
                     (8./9.)*(alp_mass*alp_mass - (4./9.)*pion_mass*pion_mass) / (alp_mass*alp_mass - eta_mass * eta_mass) /* eta */ +
                     (7./9.)*(alp_mass*alp_mass - (16./9.)*pion_mass*pion_mass) / (alp_mass*alp_mass - etap_mass*etap_mass) /* eta' */); 

  double elec_mass = Constants::Instance().elec_mass;
  double muon_mass = Constants::Instance().muon_mass;
  double tau_mass = Constants::Instance().tau_mass;
  std::complex<double> clep = 2*cAl*(B1(4*elec_mass*elec_mass / (alp_mass*alp_mass)) +
                       B1(4*muon_mass*muon_mass / (alp_mass*alp_mass)) +
                       B1(4*tau_mass*tau_mass / (alp_mass*alp_mass)));

  std::complex<double> cgamma = (5./3.)*cB + cW + cgluon + clep;

  double fsc = Constants::Instance().fine_structure_constant;
  double width = fsc*fsc*std::norm(cgamma)*alp_mass*alp_mass*alp_mass/(256*M_PI*M_PI*M_PI*fa*fa);

  return width;
}

// a -> pi+pi-pi0 
double AmplitudePiPPiMPi0(double alp_mass, double fa_eff, double m2_pip_pim) {
  double deltaI = 1./3.; // isospin violation
  double pizero_mass = Constants::Instance().pizero_mass;
  double eta_mass = Constants::Instance().eta_mass;
  double etap_mass = Constants::Instance().etap_mass;
  double fpion = Constants::Instance().fpion / sqrt(2); // fpi ~ 93MeV convention

  double M2_pi0_pi0 = pizero_mass*pizero_mass;
  double M2_eta_pi0 = - M2_pi0_pi0*sqrt(2./3.);
  double M2_etap_pi0 = -M2_pi0_pi0/sqrt(3.);
  double S_eta_pi0 = M2_eta_pi0 / (eta_mass*eta_mass - pizero_mass*pizero_mass);
  double S_etap_pi0 = M2_etap_pi0 / (etap_mass*etap_mass - pizero_mass*pizero_mass);

  double expct_a_pi = (deltaI / 2.) * (alp_mass*alp_mass) / (alp_mass*alp_mass - pizero_mass*pizero_mass);
  double expct_a_eta = (alp_mass*alp_mass/sqrt(6.) - pizero_mass*pizero_mass/sqrt(24.)) / (alp_mass*alp_mass - eta_mass*eta_mass);
  double expct_a_etap = (alp_mass*alp_mass/sqrt(12.) - pizero_mass*pizero_mass/sqrt(3.)) / (alp_mass*alp_mass - etap_mass*etap_mass);


  return (1. / (3.*fa_eff*fpion))*((3*m2_pip_pim - alp_mass*alp_mass - 2*pizero_mass*pizero_mass)*expct_a_pi\
                               - deltaI*pizero_mass*pizero_mass*(1./sqrt(3.) + sqrt(2.)*S_eta_pi0 + S_etap_pi0)*(sqrt(2)*expct_a_eta + expct_a_etap)\
                               + deltaI*(pizero_mass*pizero_mass - 2*eta_mass*eta_mass)/(pizero_mass*pizero_mass-4*eta_mass*eta_mass)*\
                                   (sqrt(3.)*pizero_mass*pizero_mass*(sqrt(2.)*S_eta_pi0 + S_etap_pi0) - 3*m2_pip_pim + alp_mass*alp_mass + 3*pizero_mass*pizero_mass));
}

double PiPPiMPi0PartialWidth(double alp_mass, double cG, double fa) {
  double piplus_mass = Constants::Instance().piplus_mass;
  double pizero_mass = Constants::Instance().pizero_mass;
  double kfactor = 2.7; // Fudge factor to match eta, eta' decay widths
  double S = 1; // symmetry factor

  if (2*piplus_mass + pizero_mass >= alp_mass) return 0.;

  // setup integration
  std::array<double, 4> masses {alp_mass, piplus_mass, piplus_mass, pizero_mass};
  std::function<double (double, double)> integrand = [alp_mass, fa, cG](double m122, double m232) -> double {return AmplitudePiPPiMPi0(alp_mass, fa/cG, m122);};
  ThreeBodyIntegrator integral(masses, integrand);

  return (kfactor / (2*S*alp_mass)) * integral.Integrate(); 
}

// a -> pi0pi0pi0
double AmplitudePi0Pi0Pi0(double alp_mass, double fa_eff) {
  double deltaI = 1./3.; // isospin violation
  double pizero_mass = Constants::Instance().pizero_mass;
  double eta_mass = Constants::Instance().eta_mass;
  double etap_mass = Constants::Instance().etap_mass;
  double fpion = Constants::Instance().fpion / sqrt(2); // fpi ~ 93MeV convention

  double M2_pi0_pi0 = pizero_mass*pizero_mass;
  double M2_eta_pi0 = - M2_pi0_pi0*sqrt(2./3.);
  double M2_etap_pi0 = -M2_pi0_pi0/sqrt(3.);
  double S_eta_pi0 = M2_eta_pi0 / (eta_mass*eta_mass - pizero_mass*pizero_mass);
  double S_etap_pi0 = M2_etap_pi0 / (etap_mass*etap_mass - pizero_mass*pizero_mass);

  double expct_a_pi = (deltaI / 2.) * (alp_mass*alp_mass) / (alp_mass*alp_mass - pizero_mass*pizero_mass);
  double expct_a_eta = (alp_mass*alp_mass/sqrt(6.) - pizero_mass*pizero_mass/sqrt(24.)) / (alp_mass*alp_mass - eta_mass*eta_mass);
  double expct_a_etap = (alp_mass*alp_mass/sqrt(12.) - pizero_mass*pizero_mass/sqrt(3.)) / (alp_mass*alp_mass - etap_mass*etap_mass);

  return (pizero_mass*pizero_mass / (fa_eff*fpion)) * (expct_a_pi - deltaI*(1./sqrt(3) + sqrt(2)*S_eta_pi0 + S_etap_pi0)*(sqrt(2.)*expct_a_eta + expct_a_etap) +\
                                                     sqrt(3)*deltaI*((pizero_mass*pizero_mass - 2*eta_mass*eta_mass) / (pizero_mass*pizero_mass - 4*eta_mass*eta_mass))*\
                                                     (sqrt(2.)*S_eta_pi0 + S_etap_pi0));
}

double Pi0Pi0Pi0PartialWidth(double alp_mass, double cG, double fa) {
  double pizero_mass = Constants::Instance().pizero_mass;
  double kfactor = 2.7; // Fudge factor to match eta, eta' decay widths
  double S = 6; // symmetry factor

  if (3*pizero_mass >= alp_mass) return 0.;

  // setup integration
  std::array<double, 4> masses {alp_mass, pizero_mass, pizero_mass, pizero_mass};
  // auto integrand = [alp_mass, cG, fa](double m122, double m232) {return pow(AmplitudePi0Pi0Pi0(alp_mass, fa/cG), 2.);};
  std::function<double (double, double)> integrand = [alp_mass, fa, cG](double m122, double m232) -> double {return AmplitudePi0Pi0Pi0(alp_mass, fa/cG);};
  ThreeBodyIntegrator integral(masses, integrand);

  return (kfactor / (2*S*alp_mass)) * integral.Integrate(); 
}

// a-> pi+pi-gamma
double PiPPiMGammaGSLIntegrand(double s, void *param);
double PiPPiMGammaIntegrand(double m2pipi, double alp_mass, double expct_a_eta, double expct_a_etap);
double GS_BW_betapi(double s);
double GS_BW_Gamma(double s, double m, double G);
double GS_BW_k(double s);
double GS_BW_h(double s);
double GS_BW_hprime(double s);
double GS_BW_d(double m);
double GS_BW_f(double s, double m, double G);
std::complex<double> RhoBreitWeitner(double s);

double PiPPiMGammaPartialWidth(double alp_mass, double cG, double fa, bool a_is_etap=false) {
  double piplus_mass = Constants::Instance().piplus_mass;
  double alpha_em = Constants::Instance().fine_structure_constant;
  double pizero_mass = Constants::Instance().pizero_mass;
  double eta_mass = Constants::Instance().eta_mass;
  double etap_mass = Constants::Instance().etap_mass;

  if (2*piplus_mass >= alp_mass) return 0.;

  double fa_eff = fa/cG;

  double m2pipi_min = 4*piplus_mass*piplus_mass;
  double m2pipi_max = alp_mass*alp_mass;

  double expct_a_eta = a_is_etap ? 0. : (alp_mass*alp_mass/sqrt(6.) - pizero_mass*pizero_mass/sqrt(24.)) / (alp_mass*alp_mass - eta_mass*eta_mass);
  double expct_a_etap = a_is_etap ? 1. :(alp_mass*alp_mass/sqrt(12.) - pizero_mass*pizero_mass/sqrt(3.)) / (alp_mass*alp_mass - etap_mass*etap_mass);

  int integrator_size = 1000;
  gsl_integration_workspace *integrator = gsl_integration_workspace_alloc(integrator_size);

  gsl_function F;
  double fparms[3] {alp_mass, expct_a_eta, expct_a_etap};
  F.function = &PiPPiMGammaGSLIntegrand; 
  F.params = fparms;
  double result, error;
  gsl_integration_qags(&F, m2pipi_min, m2pipi_max, 0., 1e-7, integrator_size, integrator, &result, &error);
 
  double integral = result;
  double prefactor = 3*alpha_em*alp_mass*alp_mass*alp_mass / (2048*pow(M_PI, 6)*fa_eff*fa_eff);

  gsl_integration_workspace_free(integrator);

  /*
  double rho_mass = Constants::Instance().rho_mass;
  double rho_width = Constants::Instance().rho_width;
  std::cout << "BETAPI: " << GS_BW_betapi(0.5) << std::endl;
  std::cout << "GAMMA: " << GS_BW_Gamma(0.5, rho_mass, rho_width) << std::endl;
  std::cout << "k: " << GS_BW_k(0.5) << std::endl;
  std::cout << "h: " << GS_BW_h(0.5) << std::endl;
  std::cout << "hprime: " << GS_BW_hprime(0.5) << std::endl;
  std::cout << "d: " << GS_BW_d(rho_mass) << std::endl;
  std::cout << "f: " << GS_BW_f(0.5, rho_mass, rho_width) << std::endl;
  std::cout << "BW: " << RhoBreitWeitner(0.5) << std::endl;
  std::cout << "PREFACTOR: " << prefactor << " INTEGRAL: " << integral << " INTEGRAND 0.5: " << PiPPiMGammaIntegrand(0.5 ,alp_mass, expct_a_eta, expct_a_etap) << std::endl;
  */

  return integral*prefactor;
}

double GS_BW_betapi(double s) {
  double pimass = Constants::Instance().piplus_mass;

  return sqrt(1 - 4.*pimass*pimass / s);
}

double GS_BW_Gamma(double s, double m, double G) {
  double beta_frac = GS_BW_betapi(s) / GS_BW_betapi(m*m);
  return G*(s/(m*m))*beta_frac*beta_frac*beta_frac;
}

double GS_BW_k(double s) {
  return (1./2.)*sqrt(s)*GS_BW_betapi(s);
}

double GS_BW_h(double s) {
  double pimass = Constants::Instance().piplus_mass;

  return (2./M_PI)*(GS_BW_k(s)/sqrt(s))*std::log((sqrt(s) + 2*GS_BW_k(s)) / (2*pimass));
}

double GS_BW_hprime(double s) {
  double pimass = Constants::Instance().piplus_mass;

  return (1./(2*M_PI*s)) + (2*pimass*pimass / (M_PI*s*s*GS_BW_betapi(s))) * std::log((1+GS_BW_betapi(s))*sqrt(s)/(2*pimass));
}

double GS_BW_d(double m) {
  double pimass = Constants::Instance().piplus_mass;
  double kval = GS_BW_k(m*m);

  return (3/M_PI)*(pimass*pimass/(kval*kval))*std::log((m+2*kval)/(2*pimass)) + m/(2*M_PI*kval) - pimass*pimass*m/(M_PI*kval*kval*kval);
}

double GS_BW_f(double s, double m, double G) {
  double kval_s = GS_BW_k(s);
  double kval_m2 = GS_BW_k(m*m);
  double hval_s = GS_BW_h(s); 
  double hval_m2 = GS_BW_h(m*m);
  double hpval = GS_BW_hprime(m*m);

  return (G*m*m/(kval_m2*kval_m2*kval_m2))*(kval_s*kval_s*(hval_s-hval_m2) + (m*m-s)*kval_m2*kval_m2*hpval);
}

std::complex<double> RhoBreitWeitner(double s) {
  double rho_mass = Constants::Instance().rho_mass;
  double rho_width = Constants::Instance().rho_width;

   return (1. + GS_BW_d(rho_mass) * rho_width / rho_mass) /\
          (rho_mass*rho_mass - s + GS_BW_f(s, rho_mass, rho_width) - 1i*rho_mass*GS_BW_Gamma(s, rho_mass, rho_width)); 
}

double PiPPiMGammaIntegrand(double m2pipi, double alp_mass, double expct_a_eta, double expct_a_etap) {
  double pimass = Constants::Instance().piplus_mass;

  double g2 = 12*M_PI;

  double arhorho = expct_a_eta/sqrt(6) + expct_a_etap/sqrt(12);
  double BW2 = std::norm(RhoBreitWeitner(m2pipi));

  return g2*g2*m2pipi*BW2*arhorho*arhorho*pow((1-m2pipi/(alp_mass*alp_mass)), 3)*pow(1-4*pimass*pimass/m2pipi, 3./2.);
}

double PiPPiMGammaGSLIntegrand(double s, void *param) {
  double *d_param = (double*)param;
  return PiPPiMGammaIntegrand(s, d_param[0], d_param[1], d_param[2]);
}

void ValidateWidths() {
  // Test out partial width calculations
  std::cout << "Expected etap -> pipigamma: 6e-5. Value: " << PiPPiMGammaPartialWidth(Constants::Instance().etap_mass, 1., Constants::Instance().fpion/sqrt(2.), true) << std::endl;

  // Checks from: https://arxiv.org/pdf/2207.08448.pdf
  // Gamma-Gamma
  std::cout << "Gamma-Gamma width check. fa=1Tev, cAl=1/100, cG=cW=cB=1\n";
  std::cout << "Ma: 0.2222, Width: 1.0972e-16. Computed: " << GammaPartialWidth(0.222, 1., 1., 1., 0.01, 1e3) << std::endl;
  std::cout << "Ma: 0.2910, Width: 1.6518e-16. Computed: " << GammaPartialWidth(0.291, 1., 1., 1., 0.01, 1e3) << std::endl;
  std::cout << "Ma: 0.4002, Width: 1.5805e-17. Computed: " << GammaPartialWidth(0.4002, 1., 1., 1., 0.01, 1e3) << std::endl;
  std::cout << "Ma: 0.5006, Width: 9.4747e-15. Computed: " << GammaPartialWidth(0.5006, 1., 1., 1., 0.01, 1e3) << std::endl;
  std::cout << "Ma: 0.6013, Width: 5.6408e-14. Computed: " << GammaPartialWidth(0.6013, 1., 1., 1., 0.01, 1e3) << std::endl;
  std::cout << "Ma: 0.7003, Width: 1.7404e-14. Computed: " << GammaPartialWidth(0.7003, 1., 1., 1., 0.01, 1e3) << std::endl;
  std::cout << "Ma: 0.8011, Width: 4.0770e-15. Computed: " << GammaPartialWidth(0.8011, 1., 1., 1., 0.01, 1e3) << std::endl;

  // 2Pi-Gamma
  std::cout << "pi+ - pi- - Gamma width check. fa=1Tev, cAl=1/100, cG=cW=cB=1\n";
  std::cout << "Ma: 0.4187, Width: 1.1543e-18. Computed: " << PiPPiMGammaPartialWidth(0.4187, 1., 1e3) << std::endl;
  std::cout << "Ma: 0.5007, Width: 2.0839e-15. Computed: " << PiPPiMGammaPartialWidth(0.5007, 1., 1e3) << std::endl;
  std::cout << "Ma: 0.6013, Width: 2.6454e-14. Computed: " << PiPPiMGammaPartialWidth(0.6013, 1., 1e3) << std::endl;
  std::cout << "Ma: 0.7003, Width: 1.8644e-14. Computed: " << PiPPiMGammaPartialWidth(0.7003, 1., 1e3) << std::endl;
  std::cout << "Ma: 0.8027, Width: 1.1448e-14. Computed: " << PiPPiMGammaPartialWidth(0.8027, 1., 1e3) << std::endl;

  // 3Pi
  std::cout << "3pi width check. fa=1Tev, cAl=1/100, cG=cW=cB=1\n";
  std::cout << "Ma: 0.4184, Width: 1.0128e-16. Computed: " << (PiPPiMPi0PartialWidth(0.4184, 1., 1e3) + Pi0Pi0Pi0PartialWidth(0.4184, 1., 1e3)) << std::endl;
  std::cout << "Ma: 0.5006, Width: 2.4837e-14. Computed: " << (PiPPiMPi0PartialWidth(0.5006, 1., 1e3) + Pi0Pi0Pi0PartialWidth(0.5006, 1., 1e3)) << std::endl;
  std::cout << "Ma: 0.6013, Width: 3.4839e-14. Computed: " << (PiPPiMPi0PartialWidth(0.6013, 1., 1e3) + Pi0Pi0Pi0PartialWidth(0.6013, 1., 1e3)) << std::endl;
  std::cout << "Ma: 0.7004, Width: 2.0601e-15. Computed: " << (PiPPiMPi0PartialWidth(0.7004, 1., 1e3) + Pi0Pi0Pi0PartialWidth(0.7004, 1., 1e3)) << std::endl;
  std::cout << "Ma: 0.8027, Width: 1.0128e-16. Computed: " << (PiPPiMPi0PartialWidth(0.8027, 1., 1e3) + Pi0Pi0Pi0PartialWidth(0.8027, 1., 1e3)) << std::endl;

}


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
  ValidateWidths();

  fReferenceRayLength = pset.get<double>("ReferenceRayLength", -1);
  fReferenceRayDistance = pset.get<double>("ReferenceRayDistance", 0.);
  fReferenceALPMass = pset.get<double>("ReferenceALPMass", -1);
  fReferenceALPcAl = pset.get<double>("ReferenceALPcAl", -1);
  fReferenceALPcG = pset.get<double>("ReferenceALPcG", -1);
  fReferenceALPcB = pset.get<double>("ReferenceALPcB", -1);
  fReferenceALPcW = pset.get<double>("ReferenceALPcW", -1);
  fReferenceALPDecayConstant = pset.get<double>("ReferenceALPDecayConstant", -1);
  fReferenceALPEnergy = pset.get<double>("ReferenceALPEnergy", -1.);
  fAllowEMDecay = pset.get<bool>("AllowEMDecay");

  // if configured to, divide out some of the decay weight
  if (fReferenceRayLength > 0. && fReferenceALPMass > 0. && fReferenceALPDecayConstant >= 0. && 
      fReferenceALPcAl >= 0. && fReferenceALPcG >= 0. && fReferenceALPcB >= 0. && fReferenceALPcW >= 0. &&
      fReferenceALPEnergy > 0.) {

    // Get each partial width
    double width_muon = MuonPartialWidth(fReferenceALPMass, fReferenceALPcAl, fReferenceALPDecayConstant);
    double width_gamma = GammaPartialWidth(fReferenceALPMass, fReferenceALPcB, fReferenceALPcW, fReferenceALPcG, fReferenceALPcAl, fReferenceALPDecayConstant);
    double width_2pig = PiPPiMGammaPartialWidth(fReferenceALPMass, fReferenceALPcG, fReferenceALPDecayConstant);
    double width_3pi = PiPPiMPi0PartialWidth(fReferenceALPMass, fReferenceALPcG, fReferenceALPDecayConstant) +\
                       Pi0Pi0Pi0PartialWidth(fReferenceALPMass, fReferenceALPcG, fReferenceALPDecayConstant);

    std::cout << "MU WIDTH: " << width_muon << " GAMMA WIDTH: " << width_gamma << std::endl;
    double partial_width = width_muon + width_gamma*fAllowEMDecay;
    double total_width = width_muon + width_gamma + width_3pi + width_2pig;

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
  double cG = flux.C3;
  double cB = flux.C4;
  double cW = flux.C5;

  // Get each partial width
  double width_muon = MuonPartialWidth(flux.mass, cAl, fA);
  double width_gamma = GammaPartialWidth(flux.mass, cB, cW, cG, cAl, fA);
  // TODO: allow decays to hadronic final states
  double width_2pig = PiPPiMGammaPartialWidth(flux.mass, cG, fA); 
  double width_3pi = PiPPiMPi0PartialWidth(flux.mass, cG, fA) + Pi0Pi0Pi0PartialWidth(flux.mass, cG, fA); 

  double total_width = width_muon + width_gamma + width_2pig + width_3pi;
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
