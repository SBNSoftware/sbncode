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
#include "sbncode/EventGenerator/MeVPrtl/Products/MeVPortalFlux.h"
#include "sbncode/EventGenerator/MeVPrtl/Products/KaonParent.h"

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
 *  @brief  Kaon2HNLFlux class definiton
 */
class Kaon2HNLFlux : public IMeVPrtlFlux
{
public:
    /**
     *  @brief  Constructor
     */
    Kaon2HNLFlux(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~Kaon2HNLFlux();

    bool MakeFlux(const simb::MCFlux &flux, MeVPortalFlux &hnl, double &weight) override;
    void configure(const fhicl::ParameterSet&) override;

    float MaxWeight() override { 
      // Weight comes from the NuMi importance weight -- max is 100 
      //
      // Compute the max weight over all the possible decay modes
      std::vector<float> max_brs(9);
      for (unsigned i = 0; i < 9; i++) {
        max_brs[i] = BranchingRatio(fM, fMagUe4, fMagUm4, i+1);
      }
      return 100. * std::max_element(max_brs.begin(), max_brs.end());
    }

private:
  // config
  float fM; //!< Mass of HNL [GeV]
  float fMagUe4;
  float fMagUm4;
  bool fKDAROnly;

  // helper functions
  double hnl_momentum(double hnl_mass, int mode);
  std::array<double, 3> rand_threebody_energy(double kaon_mass, double lep_mass, double pion_mass, double hnl_mass);
  double threebody_momentum(double kaon_mass, double lep_mass, double pion_mass, double hnl_mass, double pLambda, double pXi0);

};

Kaon2HNLFlux::Kaon2HNLFlux(fhicl::ParameterSet const &pset):
  IMeVPrtlStage("Kaon2HNLFlux"), 
  IMeVPrtlFlux(pset) 
{
    this->configure(pset);

}

//------------------------------------------------------------------------------------------------------------------------------------------

Kaon2HNLFlux::~Kaon2HNLFlux()
{
}

// static constants 
//
// kaon lifetimes
static const double kplus_lifetime = 1.238e-8; // s
static const double klong_lifetime = 5.116e-8; // s

// masses
static const double kplus_mass = 0.493677; // GeV
static const double klong_mass = 0.497611; // GeV

static const double pplus_mass = 0.139570; // GeV
static const double pzero_mass = 0.134977; // GeV

static const double muon_mass = 0.1057; // BeV
static const double elec_mass = 0.000511; // GeV

// Kaon decay parameters
static const double kaonp_mup_numu = 0.6339;
static const double kaonp_pi0_ep_nue = 0.0493;
static const double kaonp_pi0_mup_numu = 0.0330;

static const double kaonl_pim_ep_nue = 0.2020;
static const double kaonl_pip_em_nue = 0.2020;
static const double kaonl_pim_mup_numu = 0.1348;
static const double kaonl_pip_mum_numu = 0.1348;

static const double klong_e3_lambda = 0.0300;
static const double klong_e3_xi0 = -0.11;

static const double klong_m3_lambda = 0.034;
static const double klong_m3_xi0 = -0.11;

static const double kplus_e3_lambda = 0.0286;
static const double kplus_e3_xi0 = -0.35;

static const double kplus_m3_lambda = 0.033;
static const double kplus_m3_xi0 = -0.35;

double Kaon2HNLFlux::hnl_momentum(double hnl_mass, int mode) {
  switch (mode) {
    case 1 /*K0L -> nue pi- e+ */:
    case 2 /*K0L -> nuebar pi+ e-*/:
      return threebody_momentum(klong_mass, elec_mass, pplus_mass, hnl_mass, klong_e3_lambda, klong_e3_xi0);
    case 3 /* K0L -> numu pi- mu+*/:
    case 4 /*K0L -> numubar pi+ mu-*/:
      return threebody_momentum(klong_mass, muon_mass, pplus_mass, hnl_mass, klong_m3_lambda, klong_m3_xi0);
    case 6  /*K+  -> nue pi0 e+*/:
    case 9  /*K-  -> nuebar pi0 e-*/:
      return threebody_momentum(kplus_mass, elec_mass, pzero_mass, hnl_mass, kplus_e3_lambda, kplus_e3_lambda);
    case 7  /*K+  -> numu pi0 mu+*/:
    case 10 /*K-  -> numubar pi0 mu-*/:
      return threebody_momentum(kplus_mass, muon_mass, pzero_mass, hnl_mass, kplus_m3_lambda, kplus_m3_lambda);
    case 5  /*K+  -> numu mu+*/: // unreachable (2body decay)
    case 8  /*K-  -> numubar mu-*/:
      return twobody_momentum(kplus_mass, muon_mass, hnl_mass); 
    default: // unreachable
      return -1.;
  }
}

// kaon -> leption + HNL
double twobody_momentum(double kaon_mass, double lep_mass, double hnl_mass) {
  if (kaon_mass - lep_mass < hnl_mass) return -1.;

  return sqrt(kaon_mass * kaon_mass * kaon_mass * kaon_mass 
    -2 * kaon_mass * kaon_mass * lep_mass * lep_mass
    -2 * kaon_mass * kaon_mass * hnl_mass * hnl_mass
       + lep_mass * lep_mass * lep_mass * lep_mass 
       + hnl_mass * hnl_mass * hnl_mass * hnl_mass 
    -2 * lep_mass * lep_mass * hnl_mass * hnl_mass) / ( 2 * kaon_mass );
}

std::array<double, 3> Kaon2HNLFlux::rand_threebody_energy(double kaon_mass, double lep_mass, double pion_mass, double hnl_mass) {
  double sumofdaughtermass = lep_mass + pion_mass + hnl_mass;

  double r1 = GetRandom();
  double r2 = GetRandom();

  double lep_ke = std::min(r1, r2) * (kaon_mass - sumofdaughtermass);
  double pion_ke = (1 - std::max(r1, r2)) * (kaon_mass - sumofdaughtermass);
  double hnl_ke = abs(r1 - r2) * (kaon_mass - sumofdaughtermass);

  return {lep_ke + lep_mass, pion_ke + pion_mass, hnl_ke + hnl_mass};
}

// TODO: account for non-zero HNL mass
double DalitzDensity(double kaon_mass, double pLambda, double pXi0, 
                     double lep_mass, double lep_energy, 
                     double pion_mass, double pion_energy,
                     double hnl_mass, double hnl_energy) {

  double epi_max = (kaon_mass * kaon_mass + pion_mass * pion_mass - (lep_mass + hnl_mass) * (lep_mass + hnl_mass)) / (2. * kaon_mass);

  double E = epi_max - pion_energy;
  double Q2 = kaon_mass * kaon_mass + pion_mass * pion_mass - 2. * kaon_mass * pion_energy;

  double F = 1.0 + pLambda*Q2 / (pion_mass * pion_mass);
  double Fmax = 1.0;
  if (pLambda  > 0.) Fmax = 1.0 + pLambda * (1.0 + (kaon_mass * kaon_mass) / (pion_mass * pion_mass));

  double Xi = pXi0 * (1.0 + pLambda*Q2/ (pion_mass * pion_mass));

  double coeffA = kaon_mass * (2.*lep_energy*hnl_energy - kaon_mass*E) + lep_mass*lep_mass*(E /4. - hnl_energy);
  double coeffB = lep_mass * lep_mass * (hnl_energy - E/2.);
  double coeffC = lep_mass * lep_mass * E / 4.;

  double RhoMax = (Fmax*Fmax)*(kaon_mass*kaon_mass*kaon_mass/8.);

  double Rho = (F*F)*(coeffA + coeffB*Xi + coeffC*Xi*Xi);
  return Rho/RhoMax;
}

// kaon -> leption + pion + HNL
double Kaon2HNLFlux::threebody_momentum(double kaon_mass, double lep_mass, double pion_mass, double hnl_mass, double pLambda, double pXi0) {
  if (kaon_mass - lep_mass - pion_mass < hnl_mass) return -1.;

  double r, w;
  double epion, elep, ehnl;
  do {
    r = GetRandom();
    std::array<double, 3> es = rand_threebody_energy(kaon_mass, lep_mass, pion_mass, hnl_mass);
    elep = es[0];
    epion = es[1];
    ehnl = es[2];
    w = DalitzDensity(kaon_mass, pLambda, pXi0, lep_mass, elep, pion_mass, epion, hnl_mass, ehnl);
  } while (r > w);
  
  return sqrt(ehnl*ehnl - hnl_mass*hnl_mass); 
} 

double BranchingRatio(double hnl_mass, double ue4, double um4, int mode) {
  double br = 0.;
  double coupling = 0.;
  double lep_M = 0.;
  switch (mode) {
    case 1 /*K0L -> nue pi- e+ */:
      br = kaonl_pim_ep_nue;
      coupling = ue4;
      lep_M = elec_mass;
      break;
    case 2 /*K0L -> nuebar pi+ e-*/:
      br = kaonl_pip_em_nue;
      coupling = ue4;
      lep_M = elec_mass;
      break;
    case 3 /* K0L -> numu pi- mu+*/:
      br = kaonl_pim_mup_numu;
      coupling = um4;
      lep_M = muon_mass;
      break;
    case 4 /*K0L -> numubar pi+ mu-*/:
      br = kaonl_pip_mum_numu;
      coupling = um4;
      lep_M = muon_mass;
      break;
    case 5  /*K+  -> numu mu+*/:
      br = kaonp_mup_numu;
      coupling = um4;
      lep_M = muon_mass;
      break;
    case 6  /*K+  -> nue pi0 e+*/:
      br = kaonp_pi0_ep_nue;
      coupling = ue4;
      lep_M = elec_mass;
      break;
    case 7  /*K+  -> numu pi0 mu+*/:
      br = kaonp_pi0_mup_numu;
      coupling = um4;
      lep_M = muon_mass;
      break;
    case 8  /*K-  -> numubar mu-*/:
      br = kaonp_mup_numu;
      coupling = um4;
      lep_M = muon_mass;
      break;
    case 9  /*K-  -> nuebar pi0 e-*/:
      br = kaonp_pi0_ep_nue;
      coupling = ue4;
      lep_M = elec_mass;
      break;
    case 10 /*K-  -> numubar pi0 mu-*/:
      br = kaonp_pi0_mup_numu;
      coupling = um4;
      lep_M = muon_mass;
      break;
    default: // unreachable
  }
  double delta_lep_M = lep_M * lep_M / (hnl_mass * hnl_mass);
  double delta_nu_M = 0.;
  double Fm = delta_lep_M + delta_nu_M - (delta_lep_M - delta_nu_M) * (delta_lep_M - delta_nu_M);
  double lambda = 1 + delta_lep_M * delta_lep_M + delta_nu_M * delta_nu_M - 2 * (delta_lep_M + delta_nu_M + delta_lep_M * delta_nu_M);
  double sqrt_lambda = sqrt(lambda);
  double rho = Fm * sqr_lambda;
  double kinematic_factor = rho / (delta_lep_M * (1. - delta_lep_M) * (1. - delta_lep_M));

  // scale the branching ratio 
  return br * coupling * kinematic_factor;
}


//------------------------------------------------------------------------------------------------------------------------------------------
void Kaon2HNLFlux::configure(fhicl::ParameterSet const &pset)
{
  fM = pset.get<float>("M");
  fMagUm4 = pset.get<float>("MagUm4");
  fMagUe4 = pset.get<float>("MagUe4");

  fKDAROnly = pset.get<bool>("KDAROnly", false);

  std::cout << "K+ branching ratio: " << fKPBR << std::endl;
  std::cout << "K0 branching ratio: " << fKLBR << std::endl;

}

bool Kaon2HNLFlux::MakeFlux(const simb::MCFlux &flux, evgen::ldm::MeVPortalFlux &hnl, double &weight) {
  // make the kaon parent
  evgen::ldm::KaonParent kaon;
  bool success = evgen::ldm::MakeKaonParent(flux, kaon);
  if (!success) return false;

  // select on the kaon
  if (fKDAROnly && (kaon.mom.P() > 1e-3 || kaon.pos.Z() < 72000.)) return false;
  if (fKDAROnly) std::cout << "FOUND KDAR\n";

  TLorentzVector Beam4 = BeamOrigin();

  // get position in detector frame
  hnl.pos_beamcoord = kaon.pos;
  hnl.pos = kaon.pos;
  hnl.pos.Transform(fBeam2Det);
  hnl.pos += Beam4;

  // get the momentum direction in the kaon parent rest frame
  float kaon_mass = kaon.mom.M();  
  float hnl_mass = fM;
  float pion_mass = TDatabasePDG::Instance()->GetParticle(kaon.pion_pdg)->Mass();

  // ignore if we can't make this hnl
  if (kaon_mass - pion_mass < hnl_mass) return false;

  if (fRegenKinematics) {
    // new momentum!
    double p = hnl_momentum(hnl_mass, kaon.mode); 
    double e = sqrt(p*p + hnl_mass * hnl_mass);
    hnl.mom = TLorentzVector(e, p * RandomUnitVector());
  }
  else {
    // keep the Energy, rescale the momentum
    double e = flux.fnenergy;
    if (e < hnl_mass) return false;
    double p = sqrt(e * e - hnl_mass * hnl_mass); 
    TVector3 dir(flux.ndxdz * flux.fnpz / e, flux.ndydz * flux.fnpz / e, flux.fnpz / e);
    hnl.mom = TLorentzVector(e, p * dir); 
  }

  // boost to lab frame
  TLorentzVector mom;
  mom.SetVectM(kaon_frame_momentum, hnl_mass);
  mom.Boost(kaon.mom.BoostVector());

  hnl.mom_beamcoord = mom;
  // rotate to detector frame
  hnl.mom = mom;
  hnl.mom.Transform(fBeam2Det);

  hnl.kmom_beamcoord = kaon.mom;
  // also save the kaon momentum in the detector frame
  hnl.kmom = kaon.mom;
  hnl.kmom.Transform(fBeam2Det);

  // The weight is the importance weight times the branching-ratio weight 
  weight = kaon.weight * BranchingRatio(fM, fMagUe4, fMagUm4, kaon.mode);

  // set the mixing
  hnl.C1 = fMagUe4;
  hnl.C2 = fMagUm4;
  hnl.mass = fM;

  hnl.kaon_pdg = kaon.kaon_pdg;
  higgs.generator = 1; // kHNL

  return true;
}

DEFINE_ART_CLASS_TOOL(Kaon2HNLFlux)

} // namespace ldm
} // namespace evgen
