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
     *      https://arxiv.org/abs/0901.3589
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
      double MaxWeight() override { 
	return fMaxWeight; 
      }

    private:
  
      bool fVerbose;

      // Internal data
      double fMaxWeight;
      gsl_integration_workspace *fIntegrator;
      unsigned fIntegratorSize;

      // Configure the MaxWeight
      double fReferenceUE4;
      double fReferenceUM4;
      double fReferenceUT4;
      double fReferenceHNLMass;
      double fReferenceRayLength;
      double fReferenceRayDistance;
      double fReferenceHNLEnergy;
      double fReferenceHNLKaonEnergy;
      double fReferenceHNLTauEnergy;

      // Guardrail for small decay lengths
      double fMinDetectorDistance;

      // Configure the particle
      bool fMajorana;

      // Internal struct for holding decay information
      struct DecayFinalState {
	double width;
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
      typedef double(HNLMakeDecay::*HNLWidthFunction)(double hnl_mass, double ue4, double um4, double ut4);

      std::map<std::string, HNLDecayFunction> fAvailableDecays;
      std::map<std::string, HNLWidthFunction> fAvailableWidths;
      std::map<std::string, double> fAvailableDecayMasses;
      std::vector<std::string> fDecayConfig;
      std::vector<std::string> fWidthConfig;
      std::vector<HNLDecayFunction> fSelectedDecays;
      std::vector<HNLWidthFunction> fSelectedWidths;
      std::vector<HNLWidthFunction> fAllWidths;

      double TotalWidth(double hnl_mass, double ue4, double um4, double ut4);
      double SelectedWidth(double hnl_mass, double ue4, double um4, double ut4);

      double TotalWidth(const MeVPrtlFlux &flux);
      double SelectedWidth(const MeVPrtlFlux &flux);

      // Helper functions
      double CalculateMaxWeight();
      ThreebodyMomentum isotropic_threebody_momentum(double parent_mass, double childA_mass, double childB_mass, double childC_mass);
      double I1(double x, double y, double z);
      double I2(double x, double y, double z);
      double NuDiLepDecayWidth(double hnl_mass, double u4, int nu_pdg, int lep_pdg);
      double TriNuDecayWidth(double hnl_mass, double u4tot);
      double NuP0DecayWidth(double hnl_mass, double u4tot, double m0_mass, double m0_decay_const);
      double NuV0DecayWidth(double hnl_mass, double u4tot, double m0_mass, double m0_g_const);
      double LepPiWidth(double hnl_mass, double u4, double lep_mass);
      double Nul1l2Width(double hnl_mass, double ue4, double um4, double ut4,int lepplus_pdg,int lepminus_pdg);
     
	int GetWeightedNuPDG(double ue4, double um4, double ut4);

      // Width implementation functions
      double MuPiWidth(double hnl_mass, double ue4, double um4, double ut4);
      double EPiWidth(double hnl_mass, double ue4, double um4, double ut4);
      double NuMuMuWidth(double hnl_mass, double ue4, double um4, double ut4);
      double NuMuEWidth(double hnl_mass, double ue4, double um4, double ut4);
      double NuEEWidth(double hnl_mass, double ue4, double um4, double ut4);
      double TriNuWidth(double hnl_mass, double ue4, double um4, double ut4);
      double NuPi0Width(double hnl_mass, double ue4, double um4, double ut4);
      double NuEtaWidth(double hnl_mass, double ue4, double um4, double ut4);
      double NuEtaPWidth(double hnl_mass, double ue4, double um4, double ut4);
      double NuRho0Width(double hnl_mass, double ue4, double um4, double ut4);
      

      // Decay implementation functions
      DecayFinalState NuDiLep(const MeVPrtlFlux &flux, bool is_muon);
      DecayFinalState NuMupMum(const MeVPrtlFlux &flux) { return NuDiLep(flux, true); }
      DecayFinalState NuEpEm(const MeVPrtlFlux &flux) { return NuDiLep(flux, false); }
      DecayFinalState LepPi(const MeVPrtlFlux &flux, bool is_muon);
      DecayFinalState EPi(const MeVPrtlFlux &flux) { return LepPi(flux, false); }
      DecayFinalState MuPi(const MeVPrtlFlux &flux) { return LepPi(flux, true); }
      DecayFinalState NuP0(const MeVPrtlFlux &flux, int meson_pdg);
      DecayFinalState NuPi0(const MeVPrtlFlux &flux) { return NuP0(flux, 111); }
      DecayFinalState NuEta(const MeVPrtlFlux &flux) { return NuP0(flux, 221); }
      DecayFinalState NuEtaP(const MeVPrtlFlux &flux) { return NuP0(flux, 331); }
    };

    // helpers
    double lambda(double a, double b, double c) {
      return a*a + b*b + c*c - 2*a*b - 2*b*c - 2*c*a;
    }

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

    double HNLMakeDecay::MuPiWidth(double hnl_mass, double ue4, double um4, double ut4) {
      return LepPiWidth(hnl_mass, um4, Constants::Instance().muon_mass);
    }
    double HNLMakeDecay::EPiWidth(double hnl_mass, double ue4, double um4, double ut4) {
      return LepPiWidth(hnl_mass, ue4, Constants::Instance().elec_mass);
    }
    double HNLMakeDecay::NuMuMuWidth(double hnl_mass, double ue4, double um4, double ut4) {
      return NuDiLepDecayWidth(hnl_mass, ue4, 12, 13) + NuDiLepDecayWidth(hnl_mass, um4, 14, 13) + NuDiLepDecayWidth(hnl_mass, ut4, 16, 13);
    }
    double HNLMakeDecay::NuEEWidth(double hnl_mass, double ue4, double um4, double ut4) {
      return NuDiLepDecayWidth(hnl_mass, ue4, 12, 11) + NuDiLepDecayWidth(hnl_mass, um4, 14, 11) + NuDiLepDecayWidth(hnl_mass, ut4, 16, 11);
    }
    double HNLMakeDecay::TriNuWidth(double hnl_mass, double ue4, double um4, double ut4) {
      return TriNuDecayWidth(hnl_mass, ue4 + um4 + ut4);
    }
    double HNLMakeDecay::NuPi0Width(double hnl_mass, double ue4, double um4, double ut4) {
      return NuP0DecayWidth(hnl_mass, ue4 + um4 + ut4, Constants::Instance().pizero_mass, Constants::Instance().fpion); 
    }
    double HNLMakeDecay::NuEtaWidth(double hnl_mass, double ue4, double um4, double ut4) {
      return NuP0DecayWidth(hnl_mass, ue4 + um4 + ut4, Constants::Instance().eta_mass, Constants::Instance().feta); 
    }
    double HNLMakeDecay::NuEtaPWidth(double hnl_mass, double ue4, double um4, double ut4) {
      return NuP0DecayWidth(hnl_mass, ue4 + um4 + ut4, Constants::Instance().etap_mass, Constants::Instance().fetap); 
    }
    double HNLMakeDecay::NuRho0Width(double hnl_mass, double ue4, double um4, double ut4) {
      return NuV0DecayWidth(hnl_mass, ue4 + um4 + ut4, Constants::Instance().rho_mass, Constants::Instance().grho); 
    }
    double HNLMakeDecay::NuMuEWidth(double hnl_mass, double ue4, double um4, double ut4) {
      return Nul1l2Width(hnl_mass, ue4,um4,ut4,11,13)+Nul1l2Width(hnl_mass, ue4,um4,ut4,13,11); 
    }

    double HNLMakeDecay::Nul1l2Width(double hnl_mass, double ue4, double um4, double ut4,int lepplus_pdg,int lepminus_pdg) {
      double hnl_mass_pow5 = hnl_mass*hnl_mass*hnl_mass*hnl_mass*hnl_mass;
      double lepplus_mass=0;
      double lepminus_mass=0;

      if(lepminus_pdg==std::abs(11)) lepminus_mass= Constants::Instance().elec_mass;
      if(lepminus_pdg==std::abs(13)) lepminus_mass= Constants::Instance().muon_mass;

      if(lepplus_pdg==std::abs(11)) lepplus_mass= Constants::Instance().elec_mass;
      if(lepplus_pdg==std::abs(13)) lepplus_mass= Constants::Instance().muon_mass;;
      
      double u4minus=0;
      if(lepminus_pdg==std::abs(11)) u4minus=ue4;
      if(lepminus_pdg==std::abs(13)) u4minus=um4;
      if(lepminus_pdg==std::abs(15)) u4minus=ut4;

      double Gfermi = Constants::Instance().Gfermi;

      double I1val1 = I1(lepminus_mass / hnl_mass,0, lepplus_mass / hnl_mass);

      double width = (Gfermi*Gfermi*hnl_mass_pow5) * (u4minus * I1val1)/(192*M_PI*M_PI*M_PI);

      if (fMajorana) width *= 2;

      return width;


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
      double gammaC = fmod(gammaB + 2*M_PI, 2*M_PI) - M_PI;

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

      return 12.*(s - x*x - y*y)*(1 + z*z - s)*sqrt(lambda(s,x*x,y*y)*lambda(1.,s,z*z))/s;
    }

    double I2_integrand(double s, void *param) {
      double *xyz = (double *)param;
      double x = xyz[0];
      double y = xyz[1];
      double z = xyz[2];

      return 24*y*z*(1. + x*x - s)*sqrt(lambda(s,y*y,z*z)*lambda(1,s,x*x))/s;
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

    double HNLMakeDecay::TriNuDecayWidth(double hnl_mass, double u4tot) {
      double Gfermi = Constants::Instance().Gfermi;
      double hnl_mass_pow5 = hnl_mass*hnl_mass*hnl_mass*hnl_mass*hnl_mass;

      double width=Gfermi*Gfermi*hnl_mass_pow5*u4tot / (96*M_PI*M_PI*M_PI);
      //double width=Gfermi*Gfermi*hnl_mass_pow5*u4tot / (192*M_PI*M_PI*M_PI);

      return width;
    }

    // double I3(double x, double y) {
    //   return (1+2*y)*(1-y)*sqrt(lambda(1,x,y));
    // }

    double HNLMakeDecay::NuV0DecayWidth(double hnl_mass, double u4tot, double m0_mass, double m0_g_const) {
 
      if (m0_mass > hnl_mass) {
	return 0;
      } 
      
      double Gfermi = Constants::Instance().Gfermi;
      double hnl_mass_pow3 = hnl_mass*hnl_mass*hnl_mass;
      
      double mu_m0 = m0_mass*m0_mass / hnl_mass*hnl_mass;

      double width=((u4tot*Gfermi*Gfermi*hnl_mass_pow3*m0_g_const*m0_g_const) / (16*M_PI*m0_mass*m0_mass) * (1+2*m0_mass*m0_mass/(hnl_mass*hnl_mass)) * (1-mu_m0)*(1-mu_m0));

      return width; 
      
    }
    
    double HNLMakeDecay::NuDiLepDecayWidth(double hnl_mass, double u4, int nu_pdg, int lep_pdg) {
      double hnl_mass_pow5 = hnl_mass*hnl_mass*hnl_mass*hnl_mass*hnl_mass;
      double lep_mass = (lep_pdg == 13) ? Constants::Instance().muon_mass : Constants::Instance().elec_mass;

      double Gfermi = Constants::Instance().Gfermi;
      double gL = Constants::Instance().gL;
      double gR = Constants::Instance().gR;

      if (hnl_mass < lep_mass * 2.) return 0.;

      int CC = (lep_pdg+1 == nu_pdg);

      double I1val = I1(0., lep_mass / hnl_mass, lep_mass / hnl_mass);
      double I2val = I2(0., lep_mass / hnl_mass, lep_mass / hnl_mass);

      double width = (Gfermi*Gfermi*hnl_mass_pow5) * u4 * ((gL*gR/*NC*/ + CC*gR/*CC*/)*I2val + (gL*gL+gR*gR+CC*(1+2.*gL))*I1val)/(96*M_PI*M_PI*M_PI);
      //double width = (Gfermi*Gfermi*hnl_mass_pow5) * u4 * ((gL*gR/*NC*/ + CC*gR/*CC*/)*I2val + (gL*gL+gR*gR+CC*(1+2.*gL))*I1val)/(192*M_PI*M_PI*M_PI);

      return width;
    }

    HNLMakeDecay::DecayFinalState HNLMakeDecay::NuDiLep(const MeVPrtlFlux &flux, bool is_muon) {
      HNLMakeDecay::DecayFinalState ret;
      double lep_mass = is_muon ? Constants::Instance().muon_mass : Constants::Instance().elec_mass;
      int lep_pdg = is_muon ? 13 : 11;

      // Decay not kinematically allowed
      if (2*lep_mass > flux.mass) {
	ret.width = 0.;
	return ret;
      } 

      double ue4 = flux.C1;
      double um4 = flux.C2;
      double ut4 = flux.C3;
      double nue_width = NuDiLepDecayWidth(flux.mass, ue4, 12, lep_pdg);
      double numu_width = NuDiLepDecayWidth(flux.mass, um4, 14, lep_pdg);
      double nut_width = NuDiLepDecayWidth(flux.mass, ut4, 16, lep_pdg);
      double total_width = nue_width + numu_width + nut_width;

      ret.width = total_width;
      if (ret.width == 0.) return ret;

      // Three body decay
      //
      // TODO: account for anisotropies in decay
      ThreebodyMomentum momenta = isotropic_threebody_momentum(flux.mass, 0., lep_mass, lep_mass); 

      // Boost it!
      momenta.A.Boost(flux.mom.BoostVector());
      momenta.B.Boost(flux.mom.BoostVector());
      momenta.C.Boost(flux.mom.BoostVector());

      // pick whether the neutrino is nue or numu
      int nu_pdg_sign;
      if (fMajorana) {
	nu_pdg_sign = (GetRandom() > 0.5) ? 1:-1;
      }
      else {
	// same as the HNL
	nu_pdg_sign = (flux.secondary_pdg > 0) ? -1 : 1;
      }
      double nu_pdg_r = GetRandom();
      int nu_pdg = 12;
      if (nu_pdg_r > (numu_width + nue_width) / total_width) {
	nu_pdg = 16;
      }
      else if (nu_pdg_r > nue_width / total_width) {
	nu_pdg = 14;
      }
      nu_pdg = nu_pdg * nu_pdg_sign;

      ret.pdg.push_back(nu_pdg);
      ret.mom.push_back(momenta.A);

      ret.pdg.push_back(lep_pdg);
      ret.mom.push_back(momenta.B);
      ret.pdg.push_back(lep_pdg*-1);
      ret.mom.push_back(momenta.C);

      return ret;
    }

    double HNLMakeDecay::LepPiWidth(double hnl_mass, double u4, double lep_mass) {
      double piplus_mass = Constants::Instance().piplus_mass;
      double Gfermi = Constants::Instance().Gfermi;
      double fpion = Constants::Instance().fpion;
      double abs_Vud_squared = Constants::Instance().abs_Vud_squared;

      // Decay not kinematically allowed
      if (lep_mass + piplus_mass > hnl_mass) {
	return 0.;
      }

      double lep_ratio = (lep_mass * lep_mass) / (hnl_mass * hnl_mass);
      double pion_ratio = (piplus_mass * piplus_mass) / (hnl_mass * hnl_mass);
      double Ifunc = ((1+lep_ratio-pion_ratio)*(1+lep_ratio)-4*lep_ratio) * sqrt(lambda(1.,pion_ratio,lep_ratio));
      //double Ifunc = (1-pion_ratio-lep_ratio*(2+pion_ratio-lep_ratio)) * sqrt(lambda(1.,pion_ratio,lep_ratio));
      //double Ifunc = ((1-lep_ratio)*(1-lep_ratio)-pion_ratio*(1+lep_ratio)) * sqrt(lambda(1.,lep_ratio,pion_ratio));


      double width = u4 * (Gfermi * Gfermi *fpion * fpion * abs_Vud_squared * hnl_mass * hnl_mass * hnl_mass * Ifunc) / (16 * M_PI);
      // Majorana gets an extra factor b.c. it can go to pi+l- and pi-l+
      if (fMajorana) width *= 2;
  
      return width;
    }

    HNLMakeDecay::DecayFinalState HNLMakeDecay::LepPi(const MeVPrtlFlux &flux, bool is_muon) {
      HNLMakeDecay::DecayFinalState ret;
      double lep_mass = is_muon ? Constants::Instance().muon_mass : Constants::Instance().elec_mass;
      int lep_pdg = is_muon ? 13 : 11;
      double u4 = is_muon ? flux.C2 : flux.C1;

      ret.width = LepPiWidth(flux.mass, u4, lep_mass);

      // Decay not kinematically allowed
      if (ret.width == 0.) {
	return ret;
      }

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
      double piplus_mass = Constants::Instance().piplus_mass;
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
	  evgen::ldm::HNLLepPiLNCDalitz(flux.mmom, flux.sec, flux.mom, PI, LB):
	  evgen::ldm::HNLLepPiLNVDalitz(flux.mmom, flux.sec, flux.mom, PI, LB);

	assert(this_dalitz < dalitz_max);
	if (this_dalitz > dalitz_max) {
	  std::cerr << "VERY VERY BAD!!!! Incorrect dalitz max!!!\n";
	  std::cout << "VERY VERY BAD!!!! Incorrect dalitz max!!!\n";
	  std::cout << "PK: " << flux.mmom.E() << " " << flux.mmom.Px() << " " << flux.mmom.Py() << " " << flux.mmom.Pz() << std::endl;
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

      // lep
      ret.mom.push_back(LB);
      ret.pdg.push_back(lep_pdg*lep_pdg_sign);

      // pion
      ret.mom.emplace_back(PI);
      ret.pdg.push_back(211*lep_pdg_sign); // negative of lepton-charge has same-sign-PDG code

      return ret;
    }

    double HNLMakeDecay::NuP0DecayWidth(double hnl_mass, double u4tot, double m0_mass, double m0_decay_const) {
      
      if (m0_mass > hnl_mass) {
	return 0;
      } 

      double Gfermi = Constants::Instance().Gfermi;
      double hnl_mass_pow3 = hnl_mass*hnl_mass*hnl_mass;
      double mu_m0 = m0_mass*m0_mass/(hnl_mass*hnl_mass);

      double width=Gfermi*Gfermi*hnl_mass_pow3*m0_decay_const*m0_decay_const*u4tot*(1-mu_m0)*(1-mu_m0) / (64*M_PI);
      //double width=Gfermi*Gfermi*hnl_mass_pow3*m0_decay_const*m0_decay_const*u4tot*(1-mu_m0)*(1-mu_m0) / (32*M_PI);
      
      return width;
    }

    HNLMakeDecay::DecayFinalState HNLMakeDecay::NuP0(const MeVPrtlFlux &flux, int meson_pdg) {
      HNLMakeDecay::DecayFinalState ret;
      double meson_mass = 0;
      double meson_decay_constant = 0;

      switch (meson_pdg) {
      case 111:
	meson_mass = Constants::Instance().pizero_mass;
	meson_decay_constant = Constants::Instance().fpion;
	break;
      case 221:
	meson_mass = Constants::Instance().eta_mass;
	meson_decay_constant = Constants::Instance().feta;
	break;
      case 331:
	meson_mass = Constants::Instance().etap_mass;
	meson_decay_constant = Constants::Instance().fetap;
	break;
      default:
	std::cout << "Wrong pdg for NuP0 decay. Only 111, 221, 331 allowed" <<std::endl;
	exit(1);
      }

      double ue4 = flux.C1;
      double um4 = flux.C2;
      double ut4 = flux.C3;
      double total_u4 = ue4 + um4 + ut4; 

      ret.width =  NuP0DecayWidth(flux.mass, total_u4, meson_mass, meson_decay_constant);

      // Decay not kinematically allowed
      if (ret.width == 0.) {
	return ret;
      }

      // For Majorana: isotropic decay
      // For Dirac: the helicity of the initial neutrino determines the final state direction arxiv:1905.00284
      // TODO: Dirac case
  
      double p = evgen::ldm::twobody_momentum(flux.mass, 0, meson_mass);

      TLorentzVector NU;
      TLorentzVector P0;

      TVector3 dir = RandomUnitVector(); 

      NU = TLorentzVector(p*dir, sqrt(p*p + 0 * 0));
      P0 = TLorentzVector(-p*dir, sqrt(p*p + meson_mass * meson_mass));

      NU.Boost(flux.mom.BoostVector());
      P0.Boost(flux.mom.BoostVector());
 
      // Pick a neutrino pdg + sign 
      int nu_pdg_sign;
      if (fMajorana) {
	nu_pdg_sign = (GetRandom() > 0.5) ? 1:-1;
      }
      else {
	// same as the HNL
	nu_pdg_sign = (flux.secondary_pdg > 0) ? -1 : 1;
      }

      double nu_pdg = 12;
      double nu_pdg_r = GetRandom();
      if (nu_pdg_r > (ue4 + um4) / total_u4) {
	nu_pdg = 16;
      }
      else if (nu_pdg_r > ue4 / total_u4) {
	nu_pdg = 14;
      }
  
      nu_pdg = nu_pdg * nu_pdg_sign;

      // nu
      ret.mom.push_back(NU);
      ret.pdg.push_back(nu_pdg);

      // p0
      ret.mom.emplace_back(P0);
      ret.pdg.push_back(meson_pdg); 

      return ret;
    }
    double HNLMakeDecay::CalculateMaxWeight() {
      double ue4 = fReferenceUE4;
      double um4 = fReferenceUM4;
      double ut4 = fReferenceUT4;

      double hnl_mass = fReferenceHNLMass;
      double length = fReferenceRayDistance;
      double det_length = fReferenceRayLength;
      double E = fReferenceHNLEnergy;
      double P = sqrt(E*E - hnl_mass * hnl_mass);
      double hnl_gamma_beta = P/hnl_mass;


      double total_width = TotalWidth(hnl_mass, ue4, um4, ut4);

      double total_lifetime_ns = Constants::Instance().hbar / total_width;
      double total_mean_dist = total_lifetime_ns * hnl_gamma_beta * Constants::Instance().c_cm_per_ns;

      double partial_width = SelectedWidth(hnl_mass, ue4, um4, ut4);
  
      double partial_lifetime_ns = Constants::Instance().hbar / partial_width;
      double partial_mean_dist = partial_lifetime_ns * hnl_gamma_beta * Constants::Instance().c_cm_per_ns;

      if (fVerbose){
	std::cout << "Reference ue4: " << ue4 << std::endl;
	std::cout << "Reference um4: " << um4 << std::endl;
	std::cout << "Reference ut4: " << ut4 << std::endl;
	std::cout << "Reference Energy: " << E <<  " P: " << P << std::endl;
	std::cout << "REFERENCE ALL DECAY WIDTH: " << total_width << std::endl;
	std::cout << "REFERENCE ALL DECAY LENGTH: " << total_mean_dist << std::endl;
	std::cout << "REFERENCE SELECTED DECAY WIDTH: " << partial_width << std::endl;
	std::cout << "REFERENCE SELECTED DECAY LENGTH: " << partial_mean_dist << std::endl;
      }


      double weight = forcedecay_weight(total_mean_dist, length, length + det_length) * partial_width / total_width; 
      return weight;
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
  fVerbose = pset.get<bool>("Verbose", true);

  fIntegratorSize = 1000;

  fIntegrator = gsl_integration_workspace_alloc(fIntegratorSize);

  // Setup available decays
  fAvailableDecays["mu_pi"] = &HNLMakeDecay::MuPi;
  fAvailableDecayMasses["mu_pi"] = Constants::Instance().muon_mass + Constants::Instance().piplus_mass;
  
  fAvailableDecays["e_pi"] = &HNLMakeDecay::EPi;
  fAvailableDecayMasses["e_pi"] = Constants::Instance().elec_mass + Constants::Instance().piplus_mass;
  
  fAvailableDecays["nu_mu_mu"] = &HNLMakeDecay::NuMupMum;
  fAvailableDecayMasses["nu_mu_mu"] = Constants::Instance().muon_mass + Constants::Instance().muon_mass;
  
  fAvailableDecays["nu_e_e"] = &HNLMakeDecay::NuEpEm;
  fAvailableDecayMasses["nu_e_e"] = Constants::Instance().elec_mass + Constants::Instance().elec_mass;
  
  fAvailableDecays["nu_pi0"] = &HNLMakeDecay::NuPi0;
  fAvailableDecayMasses["nu_pi0"] = Constants::Instance().pizero_mass;
  
  fAvailableDecays["nu_eta"] = &HNLMakeDecay::NuEta;
  fAvailableDecayMasses["nu_eta"] = Constants::Instance().eta_mass;
  
  fAvailableDecays["nu_etap"] = &HNLMakeDecay::NuEtaP;
  fAvailableDecayMasses["nu_etap"] = Constants::Instance().etap_mass;

  // Setup available widths
  fAvailableWidths["mu_pi"] = &HNLMakeDecay::MuPiWidth;
  fAvailableWidths["e_pi"] = &HNLMakeDecay::EPiWidth;
  fAvailableWidths["nu_mu_mu"] = &HNLMakeDecay::NuMuMuWidth;
  fAvailableWidths["nu_e_e"] = &HNLMakeDecay::NuEEWidth;
  fAvailableWidths["nu_nu_nu"] = &HNLMakeDecay::TriNuWidth;
  fAvailableWidths["nu_pi0"] = &HNLMakeDecay::NuPi0Width;
  fAvailableWidths["nu_eta"] = &HNLMakeDecay::NuEtaWidth;
  fAvailableWidths["nu_etap"] = &HNLMakeDecay::NuEtaPWidth;
  fAvailableWidths["nu_rho0"] = &HNLMakeDecay::NuRho0Width;
  fAvailableWidths["nu_mu_e"] = &HNLMakeDecay::NuMuEWidth;

  // Select which ones are configued
  fDecayConfig = pset.get<std::vector<std::string>>("Decays");
  fWidthConfig = pset.get<std::vector<std::string>>("WidthDecays", fDecayConfig);

  for (const std::string &d: fDecayConfig) {
    if (fAvailableDecays.count(d)) {
      fSelectedDecays.push_back(fAvailableDecays.at(d));
      fSelectedWidths.push_back(fAvailableWidths.at(d));
      
      if (fVerbose) std::cout << "Selected Decay: " << d << std::endl;
    }
    else {
      std::cerr << "ERROR: Selected unavailable decay (" << d << ")" << std::endl;
    }
  }

  for (const std::string &d: fWidthConfig) {
    if (fAvailableWidths.count(d)) {
      fAllWidths.push_back(fAvailableWidths.at(d));
    }
    else {
      std::cerr << "ERROR: Selected unavailable decay (" << d << ")" << std::endl;
    }
  }

  fReferenceUE4 = pset.get<double>("ReferenceUE4");
  fReferenceUM4 = pset.get<double>("ReferenceUM4");
  fReferenceUT4 = pset.get<double>("ReferenceUT4");
  fReferenceHNLMass = pset.get<double>("ReferenceHNLMass");
  fReferenceRayLength = pset.get<double>("ReferenceRayLength");
  fReferenceRayDistance = pset.get<double>("ReferenceRayDistance");

  fReferenceHNLEnergy = pset.get<double>("ReferenceHNLEnergy", -1);
  fReferenceHNLKaonEnergy = pset.get<double>("ReferenceHNLEnergyFromKaonEnergy", -1.);
  fReferenceHNLTauEnergy = pset.get<double>("ReferenceHNLEnergyFromTauEnergy", -1.);
  if (fReferenceHNLEnergy < 0. && fReferenceHNLKaonEnergy > 0.) {
    double lep_mass = (fReferenceUE4 > 0) ? Constants::Instance().elec_mass : Constants::Instance().muon_mass;
    fReferenceHNLEnergy = forwardPrtlEnergy(Constants::Instance().kplus_mass, lep_mass, fReferenceHNLMass, fReferenceHNLKaonEnergy);
  }
  else if (fReferenceHNLEnergy < 0. && fReferenceHNLTauEnergy > 0.) {
    fReferenceHNLEnergy = forwardPrtlEnergy(Constants::Instance().tau_mass, Constants::Instance().piplus_mass, fReferenceHNLMass, fReferenceHNLTauEnergy);
  }

  fMinDetectorDistance = pset.get<double>("MinDetectorDistance", 100e2); // 100m for NuMI -> SBN/ICARUS

  fMajorana = pset.get<bool>("Majorana");

  fMaxWeight = CalculateMaxWeight();

}


double HNLMakeDecay::TotalWidth(const MeVPrtlFlux &flux) {
  return TotalWidth(flux.mass, flux.C1, flux.C2, flux.C3);
}

double HNLMakeDecay::SelectedWidth(const MeVPrtlFlux &flux) {
  return SelectedWidth(flux.mass, flux.C1, flux.C2, flux.C3);
}

double HNLMakeDecay::TotalWidth(double hnl_mass, double ue4, double um4, double ut4) {
  double ret = 0.;
  for (const HNLMakeDecay::HNLWidthFunction F: fAllWidths) {
    ret += (*this.*F)(hnl_mass, ue4, um4, ut4);
  }

  return ret;

}

double HNLMakeDecay::SelectedWidth(double hnl_mass, double ue4, double um4, double ut4) {
  double ret = 0.;
  for (const HNLMakeDecay::HNLWidthFunction F: fSelectedWidths) {
    ret += (*this.*F)(hnl_mass, ue4, um4, ut4);
  }

  return ret;

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
  double partial_width = 0.;
  for (const HNLMakeDecay::HNLDecayFunction F: fSelectedDecays) {
    decays.push_back((*this.*F)(flux)); 
    partial_width += decays.back().width;
  }
  

  if (partial_width == 0.) return false;

  // pick one
  double sum_width = 0.;
  int idecay = decays.size()-1;
  double rand = GetRandom();
  for (unsigned i = 0; i < decays.size()-1; i++) {
    sum_width += decays[i].width;
    if (rand < sum_width / partial_width) {
      idecay = i;
      break;
    }
  }

  // Get the decay probability

  // partial lifetime
  double partial_lifetime_ns = Constants::Instance().hbar / partial_width;

  // multiply by gamma*v to get the length
  double partial_mean_dist = partial_lifetime_ns * flux.mom.Gamma() * flux.mom.Beta() * Constants::Instance().c_cm_per_ns;

  double in_dist = (flux.pos.Vect() - in).Mag();
  double out_dist = (flux.pos.Vect() - out).Mag();

  // Total width
  double total_width = TotalWidth(flux);
  double total_lifetime_ns = Constants::Instance().hbar / total_width;
  double total_mean_dist = total_lifetime_ns * flux.mom.Gamma() * flux.mom.Beta() * Constants::Instance().c_cm_per_ns;

  if (fVerbose){  
  std::cout <<"Trinu Branching Ratio: " << HNLMakeDecay::TriNuWidth(flux.mass, flux.C1, flux.C2, flux.C3)/total_width << std::endl;
  std::cout <<"NuPi0 Branching Ratio: " << HNLMakeDecay::NuPi0Width(flux.mass, flux.C1, flux.C2, flux.C3)/total_width << std::endl;
  std::cout <<"mupi Branching Ratio: " << HNLMakeDecay::MuPiWidth(flux.mass, flux.C1, flux.C2, flux.C3)/total_width << std::endl;
  std::cout <<"epi Branching Ratio: " << HNLMakeDecay::EPiWidth(flux.mass, flux.C1, flux.C2, flux.C3)/total_width << std::endl;
  std::cout <<"nuMuMu Branching Ratio: " << HNLMakeDecay::NuMuMuWidth(flux.mass, flux.C1, flux.C2, flux.C3)/total_width << std::endl;
  std::cout <<"NuMuE Branching Ratio: " << HNLMakeDecay::NuMuEWidth(flux.mass, flux.C1, flux.C2, flux.C3)/total_width << std::endl;
  std::cout <<"nuEE Branching Ratio: " << HNLMakeDecay::NuEEWidth(flux.mass, flux.C1, flux.C2, flux.C3)/total_width << std::endl;
  std::cout <<"nueta Branching Ratio: " << HNLMakeDecay::NuEtaWidth(flux.mass, flux.C1, flux.C2, flux.C3)/total_width << std::endl;
  std::cout <<"nuetaP Branching Ratio: " << HNLMakeDecay::NuEtaPWidth(flux.mass, flux.C1, flux.C2, flux.C3)/total_width << std::endl;
  std::cout <<"nurho0 Branching Ratio: " << HNLMakeDecay::NuRho0Width(flux.mass, flux.C1, flux.C2, flux.C3)/total_width << std::endl;
  
  std::cout << "total Branching Ratio: " << total_width/total_width << std::endl;
  }

if (fVerbose){
    std::cout << "TOTAL DECAY WIDTH: " << total_width << std::endl;
    std::cout << "TOTAL DECAY DIST: " << total_mean_dist << std::endl;
    std::cout << "SELECTED DECAY WIDTH: " << partial_width << std::endl;
    std::cout << "SELECTED DECAY DIST: " << partial_mean_dist << std::endl;
    std::cout << "DISTANCE TO DETECTOR: " << in_dist << " DISTANCE OUT OF DETECTOR: " << out_dist << std::endl;
    std::cout << "force decay weight: " << forcedecay_weight(total_mean_dist, in_dist, out_dist) << std::endl;
    std::cout << "partial/total width: " << partial_width/total_width << std::endl;
    std::cout << std::endl;
  }  

  
   // saves the weight
  weight = forcedecay_weight(total_mean_dist, in_dist, out_dist) * partial_width / total_width; 

  // ignore events that will never reach the detector
  if (weight == 0.) return false;

  // Get the decay location 
  double flat_rand = CLHEP::RandFlat::shoot(fEngine, 0, 1.);

  double decay_rand = flat_to_exp_rand(flat_rand, total_mean_dist, in_dist, out_dist);
  TVector3 decay_pos = flux.pos.Vect() + decay_rand * (in - flux.pos.Vect()).Unit();

  // Save the decay info
  decay.pos = TLorentzVector(decay_pos, TimeOfFlight(flux, decay_pos));
  for (const TLorentzVector &p: decays[idecay].mom) {
    decay.daughter_mom.push_back(p.Vect());
    decay.daughter_e.push_back(p.E());
  }
  decay.daughter_pdg = decays[idecay].pdg;

  decay.total_decay_width = total_width;
  decay.total_mean_lifetime = total_lifetime_ns;
  decay.total_mean_distance = total_mean_dist;
  decay.allowed_decay_fraction = partial_width / total_width;

  return true;
}

DEFINE_ART_CLASS_TOOL(HNLMakeDecay)

} // namespace ldm
} // namespace evgen
