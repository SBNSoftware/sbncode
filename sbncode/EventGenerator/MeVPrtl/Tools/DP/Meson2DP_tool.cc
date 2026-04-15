

#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "sbncode/EventGenerator/MeVPrtl/Tools/IMeVPrtlFlux.h"
#include "sbnobj/Common/EventGen/MeVPrtl/MeVPrtlFlux.h"
#include "sbnobj/Common/EventGen/MeVPrtl/MesonParent.h"

#include "TDatabasePDG.h"
#include "TParticlePDG.h"

#include <cmath>
#include <algorithm>
#include <string>
#include <vector>

namespace evgen {
namespace ldm {

class Meson2DP : public IMeVPrtlFlux {
public:
  Meson2DP(fhicl::ParameterSet const& pset);
  ~Meson2DP() override = default;

  void configure(fhicl::ParameterSet const& pset) override;
  bool MakeFlux(const simb::MCFlux& flux, MeVPrtlFlux& dp, double& weight) override;

  double MaxWeight() override { return 1.0; }

private:
  double fM      = 0.0;
  double fe      = 0.0;

  bool fVerbose = false;

  std::vector<int> fParents;
  bool fCorrectDK2NU = true;

  std::string fSecondarySignConvention = "legacy"; // legacy/physical

  static double pdg_mass_GeV(int pdg) {
    auto* db = TDatabasePDG::Instance();
    TParticlePDG* p = db->GetParticle(pdg);
    return p ? p->Mass() : -1.0;
  }

  static double BRMeson2SM(int pdg){

    if (pdg == 221){
      return 0.39;
    }
    else if (pdg == 111){
      return 0.99;
    }
    return 0.0;
  }

  static double BRMeson2DP(int pdg, double e, double mA){
    
    double m_meson = pdg_mass_GeV(pdg);
    double temp1 = 2*e*e;
    double temp2 = 1-((mA*mA)/(m_meson*m_meson));
    double br = BRMeson2SM(pdg);
    double br_dp = temp1*temp2*br;

    return br_dp;

  }

  static double twobody_momentum(double M, double m1, double m2) {
    if (M < m1 + m2) return -1.0;
    const double a = M*M - (m1+m2)*(m1+m2);
    const double b = M*M - (m1-m2)*(m1-m2);
    if (a < 0 || b < 0) return -1.0;
    return std::sqrt(a*b) / (2.0*M);
  }


  bool parent_allowed(int abs_pdg) const {
    for (int p : fParents) if (abs_pdg == std::abs(p)) return true;
    return false;
  }
};

Meson2DP::Meson2DP(fhicl::ParameterSet const& pset)
  : IMeVPrtlStage("Meson2DP")
  , IMeVPrtlFlux(pset)
{
  configure(pset);
}

void Meson2DP::configure(fhicl::ParameterSet const& pset)
{
  fM      = pset.get<double>("M");
  fe      = pset.get<double>("e");

  fVerbose = pset.get<bool>("Verbose", false);

  fParents = pset.get<std::vector<int>>("Parents", std::vector<int>{321}); // default K±
  fCorrectDK2NU = pset.get<bool>("CorrectDK2NU", true);

  fSecondarySignConvention =
    pset.get<std::string>("SecondarySignConvention", "legacy");
  std::transform(fSecondarySignConvention.begin(), fSecondarySignConvention.end(),
                 fSecondarySignConvention.begin(), ::tolower);
}

bool Meson2DP::MakeFlux(const simb::MCFlux& flux,
                            evgen::ldm::MeVPrtlFlux& dp,
                            double& weight)
{
  evgen::ldm::MesonParent parent(flux);

  const int mpdg = parent.meson_pdg;
  const int abs_mpdg = std::abs(mpdg);

  if (!parent_allowed(abs_mpdg)) return false;


  if (abs_mpdg != 221) return false; // eta only

  const double mA = fM;
  if (!(mA > 0.0)) return false;

  const double e = fe;
  if (!(e > 0.0)) return false;


  double br = BRMeson2DP(abs_mpdg, e, mA);
  if (!(br > 0.0)) return false;
  if (br > 1.0) {
    if (fVerbose) mf::LogWarning("Meson2DP") << "BRprod>1 (" << br << "), capping to 1.";
    br = 1.0;
  }

  
  TLorentzVector Beam4 = BeamOrigin();
  dp.pos_beamcoord = parent.pos;
  dp.pos = parent.pos;
  dp.pos.Transform(fBeam2Det);
  dp.pos += Beam4;

  dp.mmom_beamcoord = parent.mom;
  dp.mmom = parent.mom;
  dp.mmom.Transform(fBeam2Det);


  const int gamma_pdg = 22;

  const double mMeson = pdg_mass_GeV(221);
  const double mgamma = 0.0;

  const double pstar = twobody_momentum(mMeson, mgamma, mA);
  if (pstar < 0) return false;

  const double eA = std::sqrt(pstar*pstar + mA*mA);

  // isotropic S in K rest frame
  dp.mom = TLorentzVector(pstar * RandomUnitVector(), eA);

  // boost to lab
  TLorentzVector mom_lab = dp.mom;
  mom_lab.Boost(parent.mom.BoostVector());

  dp.mom_beamcoord = mom_lab;
  dp.mom = mom_lab;
  dp.mom.Transform(fBeam2Det);

  dp.sec_beamcoord = dp.mmom_beamcoord - dp.mom_beamcoord;
  dp.sec = dp.mmom - dp.mom;

  double denom = 1.0;

  weight = parent.weight * br / denom;
  if (weight == 0.0) return false;

  dp.mass = mA;
  dp.meson_pdg = mpdg;
  dp.secondary_pdg = gamma_pdg;

  dp.generator = 1;
  dp.equiv_enu = EnuLab(flux.fnecm, dp.mmom, dp.pos);

 
  dp.C1 = 0.0; dp.C2 = 0.0; dp.C3 = 0.0; dp.C4 = 0.0; dp.C5 = 0.0;

  if (fVerbose) {
    mf::LogInfo("Meson2DP")
      << "parent=" << mpdg
      << " mA=" << mA
      << " BRprod=" << br
      << " denom=" << denom
      << " weight=" << weight;
  }

  return true;
}

DEFINE_ART_CLASS_TOOL(Meson2DP)

} // namespace ldm
} // namespace evgen
