
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "sbncode/EventGenerator/MeVPrtl/Tools/IMeVPrtlDecay.h"
#include "sbncode/EventGenerator/MeVPrtl/Tools/Constants.h"
#include "sbnobj/Common/EventGen/MeVPrtl/MeVPrtlFlux.h"

#include "CLHEP/Random/RandFlat.h"

#include <cmath>
#include <algorithm>
#include <string>
#include <fstream>
#include <sstream>

namespace evgen {
namespace ldm {

namespace {

double flat_to_exp_rand(double x, double mean, double a, double b) {
  const double A = (1. - std::exp(-(b-a)/mean));
  return - mean * std::log(1 - x * A) + a;
}

double forcedecay_weight(double mean, double a, double b) {
  return std::exp(-a/mean) - std::exp(-b/mean);
}


} // anon

class DPDecay : public IMeVPrtlDecay {
public:
  DPDecay(fhicl::ParameterSet const& pset);
  ~DPDecay() override = default;

  void configure(fhicl::ParameterSet const& pset) override;

  bool Decay(const MeVPrtlFlux& flux,
             const TVector3& in,
             const TVector3& out,
             MeVPrtlDecay& decay,
             double& weight) override;


  double MaxWeight() override { return fMaxWeight; }

private:
  bool fAllowElectronDecay = true;
  bool fAllowMuonDecay     = true;
  bool fAddTimeOfFlight    = true;
  bool fVerbose            = false;

  // FHiCL inputs
  //  double fCTau_cm = 1.0;
  //  double fBR_ee   = 1.0;
  double fEpsilon = -1.0;
  std::vector<double> fMassGrid;
  std::vector<double> fBRHad;
  std::vector<double> fBRMu;
  std::vector<double> fBRE;

  double fReferenceRayLength   = -1;
  double fReferenceRayDistance = 0.;
  double fReferenceMass        = -1;
  double fReferenceCTau_cm     = -1;
  double fReferenceEnergy      = -1;

  double fMaxWeight = -1.;
  double interp(const std::vector<double>& x,
		const std::vector<double>& y,
		double val) const;

  double Gamma_ll(double m, double eps, double ml) const;

  double Gamma_total(double m, double eps) const;

};
  
  auto load_table = [](const std::string& fname,
		       std::vector<double>& m,
		       std::vector<double>& br) {
    std::ifstream f(fname);
    double x, y;
    while (f >> x >> y) {
      m.push_back(x);
      br.push_back(y);
    }
  };


DPDecay::DPDecay(fhicl::ParameterSet const& pset)
  : IMeVPrtlStage("DPDecay")
{
  configure(pset);
}

void DPDecay::configure(fhicl::ParameterSet const& pset)
{
  fAllowElectronDecay = pset.get<bool>("AllowElectronDecay", true);
  fAllowMuonDecay     = pset.get<bool>("AllowMuonDecay", true);
  fAddTimeOfFlight    = pset.get<bool>("AddTimeOfFlight", true);
  fVerbose            = pset.get<bool>("Verbose", false);
  fEpsilon = pset.get<double>("Epsilon");
  //fCTau_cm = pset.get<double>("CTau_cm");
  //fBR_ee   = pset.get<double>("BR_ee");

  /*
  if (fCTau_cm <= 0.0) {
    throw cet::exception("DPDecay") << "CTau_cm must be > 0.";
  }
  fBR_ee = std::clamp(fBR_ee, 0.0, 1.0);
*/
  fReferenceRayLength   = pset.get<double>("ReferenceRayLength", -1.0);
  fReferenceRayDistance = pset.get<double>("ReferenceRayDistance", 0.0);
  fReferenceMass        = pset.get<double>("ReferenceMass", -1.0);
  fReferenceCTau_cm     = pset.get<double>("ReferenceCTau_cm", -1.0);
  fReferenceEnergy      = pset.get<double>("ReferenceEnergy", -1.0);


  load_table("bfrac_dark_photon_e_e.txt", fMassGrid, fBRE);
  load_table("bfrac_dark_photon_mu_mu.txt", fMassGrid, fBRMu);
  load_table("bfrac_dark_photon_hadrons.txt", fMassGrid, fBRHad);

  if (fReferenceRayLength > 0 && fReferenceMass > 0 && fReferenceCTau_cm > 0 && fReferenceEnergy > 0) {
    const double mS = fReferenceMass;

    const double p = std::sqrt(std::max(0.0, fReferenceEnergy*fReferenceEnergy - mS*mS));
    const double gamma_beta = (mS > 0) ? (p / mS) : 0.0;

    const double mean_dist = fReferenceCTau_cm * gamma_beta;
    if (mean_dist <= 0) { fMaxWeight = 0.0; return; }

    fMaxWeight = forcedecay_weight(mean_dist, fReferenceRayDistance,
                                   fReferenceRayDistance + fReferenceRayLength);
  } else {
    fMaxWeight = -1.;
  }
}

  double DPDecay::interp(const std::vector<double>& x,
		const std::vector<double>& y,
		double val) const {
    int idx = 0;
    double best = 1e9;
    for (size_t i = 0; i < x.size(); i++) {
      double d = std::abs(x[i] - val);
      if (d < best) { best = d; idx = i; }
    }
    return y[idx];
  }
  double DPDecay::Gamma_ll(double m, double eps, double ml) const {
    const double alpha = 1.0/137.0;
    if (m <= 2*ml) return 0.0;

    double beta = std::sqrt(1 - 4*ml*ml/(m*m));
    return (1.0/3.0) * alpha * eps*eps * m *
      (1 + 2*ml*ml/(m*m)) * beta;
  }
  double DPDecay::Gamma_total(double m, double eps) const {
    if (m < 0.29) {
      return Gamma_ll(m, eps, 0.000511) +
	Gamma_ll(m, eps, 0.105);
    } else {
      double br_mu = interp(fMassGrid, fBRMu, m);
      double gmu = Gamma_ll(m, eps, 0.105);
      return (br_mu > 0) ? gmu / br_mu : 0.0;
    }
  }

bool DPDecay::Decay(const MeVPrtlFlux& flux,
                        const TVector3& in,
                        const TVector3& out,
                        MeVPrtlDecay& decay,
                        double& weight)
{
  const auto& C = Constants::Instance();

  const double mS = flux.mass;
  if (!(mS > 0.0)) return false;

  if (mS <= 2.0 * C.elec_mass) {
    throw cet::exception("DPDecay")
      << "BAD MASS: mS=" << mS << " below 2*me threshold.";
  }

  /*
  const bool ee_open  = fAllowElectronDecay && (mS > 2.0*C.elec_mass);
  const bool mumu_open = fAllowMuonDecay    && (mS > 2.0*C.muon_mass);

  if (!ee_open && !mumu_open) return false;

  double br_e = 0.0;
  double br_mu = 0.0;

  if (ee_open && mumu_open) {
    br_e  = fBR_ee;
    br_mu = 1.0 - fBR_ee;
  } else if (ee_open) {
    // mu closed -> force ee
    br_e = 1.0;
    br_mu = 0.0;
  } else {
    // ee closed (shouldn't happen if mS>2me) -> force mu
    br_e = 0.0;
    br_mu = 1.0;
  }
  */
  // mean lab decay distance [cm]
  //  const double mean_dist = fCTau_cm * flux.mom.Gamma() * flux.mom.Beta();

  const double mA = flux.mass;
  const double Gamma = Gamma_total(mA, fEpsilon);

  double Gamma_ee = Gamma_ll(mA, fEpsilon, 0.000511);
  double BR_ee = Gamma_ee / Gamma;

  if (Gamma <= 0) return false;

  // hbar*c in GeV*cm
  const double hbarc = 1.973e-14;

  double ctau_cm = hbarc / Gamma;

  const double mean_dist = ctau_cm * flux.mom.Gamma() * flux.mom.Beta();
  if (mean_dist <= 0.0) return false;

  const double in_dist  = (flux.pos.Vect() - in ).Mag();
  const double out_dist = (flux.pos.Vect() - out).Mag();

  weight = forcedecay_weight(mean_dist, in_dist, out_dist);
  weight *= BR_ee;
  if (weight == 0.0) return false;

  
  // pick channel by BR
  /*  const double u = CLHEP::RandFlat::shoot(fEngine, 0.0, 1.0);
  const int ch = pick_by_br(br_e, br_mu, u);
  if (ch < 0) return false;
  */
  const int abs_pdg = 11;
  const double ml   = C.elec_mass;

  // sample decay point along ray within [in,out]
  const double u2 = CLHEP::RandFlat::shoot(fEngine, 0.0, 1.0);
  const double decay_rand = flat_to_exp_rand(u2, mean_dist, in_dist, out_dist);
  TVector3 decay_pos3 = flux.pos.Vect() + decay_rand * (in - flux.pos.Vect()).Unit();

  const double decay_time = fAddTimeOfFlight ? TimeOfFlight(flux, decay_pos3) : flux.pos.T();
  TLorentzVector decay_pos(decay_pos3, decay_time);

  // 2-body kinematics in S rest frame
  const double E_rf = mS / 2.0;
  const double p_rf = std::sqrt(std::max(0.0, E_rf*E_rf - ml*ml));

  TVector3 pA_rf = RandomUnitVector() * p_rf;
  TVector3 pB_rf = -pA_rf;

  TLorentzVector p4A(pA_rf, E_rf);
  TLorentzVector p4B(pB_rf, E_rf);

  p4A.Boost(flux.mom.BoostVector());
  p4B.Boost(flux.mom.BoostVector());


  //const double tau_ns = fCTau_cm / C.c_cm_per_ns;
  double tau_ns = (ctau_cm > 0) ? (ctau_cm / C.c_cm_per_ns) : 0.0;
  const double gamma_tot = (tau_ns > 0) ? (C.hbar / tau_ns) : 0.0;

  decay.total_decay_width   = gamma_tot;
  decay.total_mean_lifetime = tau_ns;
  decay.total_mean_distance = mean_dist;
  decay.allowed_decay_fraction = 1.0;

  decay.pos = decay_pos;

  decay.daughter_mom.clear();
  decay.daughter_e.clear();
  decay.daughter_pdg.clear();

  decay.daughter_mom.push_back(p4A.Vect());
  decay.daughter_e.push_back(p4A.E());
  decay.daughter_pdg.push_back(+abs_pdg);

  decay.daughter_mom.push_back(p4B.Vect());
  decay.daughter_e.push_back(p4B.E());
  decay.daughter_pdg.push_back(-abs_pdg);

  if (fVerbose) {
    mf::LogInfo("DPDecay")
      << "mS=" << mS
      << " BR_ee(eff)=" << BR_ee 
      << " mean_dist_cm=" << mean_dist
      << " weight=" << weight;
  }

  return true;
}

DEFINE_ART_CLASS_TOOL(DPDecay)

} // namespace ldm
} // namespace evgen
