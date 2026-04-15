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
 *  @brief  PiZero2DPFlux class definiton
 *
 *  Implementation of Kaon->DP branching ratio taken from:
 *      arXiv:1912.07622
 */
class PiZero2DPFlux : public IMeVPrtlFlux
{
public:
    /**
     *  @brief  Constructor
     */
    PiZero2DPFlux(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~PiZero2DPFlux();

    bool MakeFlux(const simb::MCFlux &flux, MeVPrtlFlux &dp, double &weight) override;
    void configure(const fhicl::ParameterSet&) override;

    double MaxWeight() override;

private:
  // config
  double fM; //!< Mass of DP [GeV]
  double fe;
  double fTarget2Absorber;
  double Pi0ToDPBR(double mA) const;
};

PiZero2DPFlux::PiZero2DPFlux(fhicl::ParameterSet const &pset):
  IMeVPrtlStage("PiZero2DPFlux"), 
  IMeVPrtlFlux(pset) 
{
    this->configure(pset);

}

//------------------------------------------------------------------------------------------------------------------------------------------

PiZero2DPFlux::~PiZero2DPFlux()
{
}

// helper functions
//
// pi0 -> photon + DP
double dp_momentum(double pi0_mass, double gamma_mass, double dp_mass) {
  if (pi0_mass - gamma_mass < dp_mass) return -1.;

  return sqrt(pi0_mass * pi0_mass * pi0_mass * pi0_mass 
    -2 * pi0_mass * pi0_mass * gamma_mass * gamma_mass
    -2 * pi0_mass * pi0_mass * dp_mass * dp_mass
       + gamma_mass * gamma_mass * gamma_mass * gamma_mass 
       + dp_mass * dp_mass * dp_mass * dp_mass 
    -2 * gamma_mass * gamma_mass * dp_mass * dp_mass) / ( 2 * pi0_mass );
}


double SMKaonBR(int pi0_pdg) {
  // The Kaons in Dk2nu file only include those that decay to neutrinos.
  //
  // We want all kaons -- in order to account for this, we divide by the 
  // branching-ratio of kaons to neutrinos
  //
  // Taken from: 
  // /cvmfs/minerva.opensciencegrid.org/minerva/beamsim/x86_64/geant4/source/particles/hadrons/mesons/src/G4KaonPlus.cc
  // /cvmfs/minerva.opensciencegrid.org/minerva/beamsim/x86_64/geant4/source/particles/hadrons/mesons/src/G4KaonZerLong.cc
  switch (pi0_pdg) {
    case 321:
      return 0.6339 /* 5 */ + 0.0559 /* 6 */ + 0.0330 /* 7 */;
    case -321:
      return 0.6339 /* 8 */ + 0.0559 /* 9 */ + 0.0330 /* 10 */;
    case 130:
      return 0.2020 /* 1 */ + 0.2020 /* 2 */ + 0.1348 /* 3 */ + 0.1348 /* 4 */;
    default: 
      return -1;
  }
}

double BranchingRatio(double dp_mass, double u4, bool is_muon) {
  double lep_mass = is_muon ? Constants::Instance().muon_mass : Constants::Instance().elec_mass;
  double kplus_mass = Constants::Instance().kplus_mass;

  if (dp_mass > kplus_mass - lep_mass) return 0.;

  double smbr = is_muon ? Constants::Instance().kaonp_mup_numu : Constants::Instance().kaonp_ep_nue;
  double lep_ratio = (lep_mass * lep_mass) / (kplus_mass * kplus_mass);
  double dp_ratio = (dp_mass * dp_mass) / (kplus_mass * kplus_mass);
  double kinematic_factor = (lep_ratio + dp_ratio - (lep_ratio - dp_ratio) * (lep_ratio - dp_ratio)) \
       * sqrt(1 + dp_ratio * dp_ratio +  lep_ratio * lep_ratio - 2*(dp_ratio + lep_ratio + dp_ratio*lep_ratio)) \
       / (lep_ratio * (1. - lep_ratio) * (1. - lep_ratio));

  // scale the branching ratio 
  return smbr * (u4 / (1. - u4)) * kinematic_factor;
}

std::pair<double, bool> Branch(double dp_mass, double ue4, double um4, double rand) {
  double kplus_mass = Constants::Instance().kplus_mass;

  double br_muon = (um4 > 0. && dp_mass < kplus_mass - Constants::Instance().muon_mass) ? BranchingRatio(dp_mass, um4, true) : 0.;
  double br_elec = (ue4 > 0. && dp_mass < kplus_mass - Constants::Instance().elec_mass) ? BranchingRatio(dp_mass, ue4, false): 0.;
  if (br_muon == 0. && br_elec == 0.) return {0., false};

  double br_weight = br_muon + br_elec;

  bool is_muon = rand < (br_muon / br_weight);

  return {br_weight, is_muon};
}

  double PiZero2DPFlux::Pi0ToDPBR(double mA) const {

    double mpi0 = Constants::Instance().pizero_mass;
    double br_gammagamma = 0.988; // PDG value

    if (mA >= mpi0) return 0.0;

    double factor = 1.0 - (mA * mA) / (mpi0 * mpi0);
    return 2.0 * fe * fe * std::pow(factor, 3) * br_gammagamma;
  }

//------------------------------------------------------------------------------------------------------------------------------------------
void PiZero2DPFlux::configure(fhicl::ParameterSet const &pset)
{
  fM = pset.get<double>("M");
  fe = pset.get<double>("e");
  //  fMagUm4 = pset.get<double>("MagUm4");
  //  fMagUe4 = pset.get<double>("MagUe4");

  fTarget2Absorber = pset.get<double>("Target2Absorber", 5000);
  //  fKDAROnly = pset.get<bool>("KDAROnly", false);
  /*
  double max_mass = (fMagUe4 > 0.) ? (Constants::Instance().kplus_mass - Constants::Instance().elec_mass) : 
      (Constants::Instance().kplus_mass - Constants::Instance().muon_mass);

  if (fM > max_mass) {
    throw cet::exception("PiZero2DPFlux Tool: BAD MASS. Configured mass (" + std::to_string(fM) +
         ") is larger than maximum allowed by enabled couplings (" + std::to_string(max_mass) +  ").");
  }
  */
}

double PiZero2DPFlux::MaxWeight() { 
  // Weight comes from the NuMi importance weight -- max is 100 (add in an epsilon)
  // Scale by the branching ratios here
  //  return 100.0001 * std::max(BranchingRatio(fM, fMagUe4, false) / SMKaonBR(321), BranchingRatio(fM, fMagUm4, true) / SMKaonBR(321));
  return 1;
}

bool PiZero2DPFlux::MakeFlux(const simb::MCFlux &flux, evgen::ldm::MeVPrtlFlux &dp, double &weight) {
  // make the pi0 parent
  evgen::ldm::MesonParent kaon(flux);
  if (abs(kaon.meson_pdg) != 111) return false; // Only take charged kaons
  double mpi0 = Constants::Instance().pizero_mass;
  TLorentzVector Beam4 = BeamOrigin();

  // get position in detector frame
  dp.pos_beamcoord = kaon.pos;
  dp.pos = kaon.pos;
  dp.pos.Transform(fBeam2Det);
  dp.pos += Beam4;

  // Branch the parent Kaon Decay
  double dp_mass = fM;
  //  std::pair<double, bool> decay = Branch(dp_mass, fMagUe4, fMagUm4, GetRandom());
  //double br = decay.first;

  double br = Pi0ToDPBR(fM);
  
  //bool is_muon = decay.second;
  //double lep_mass = is_muon ? Constants::Instance().muon_mass : Constants::Instance().elec_mass;
  if (br == 0.) return false;
  if (fVerbose) std::cout << "BR: " << br << std::endl;

  // get the momentum direction in the kaon parent rest frame
  double p = (mpi0 * mpi0 - fM * fM) / (2.0 * mpi0);
  double e = sqrt(p*p + dp_mass * dp_mass);
  // Two-body decays are isotropic
  dp.mom = TLorentzVector(p*RandomUnitVector(), e);

  // boost to lab frame
  TLorentzVector mom = dp.mom;
  mom.Boost(kaon.mom.BoostVector());

  dp.mom_beamcoord = mom;
  // rotate to detector frame
  dp.mom = mom;
  dp.mom.Transform(fBeam2Det);

  dp.mmom_beamcoord = kaon.mom;
  // also save the kaon momentum in the detector frame
  dp.mmom = kaon.mom;
  dp.mmom.Transform(fBeam2Det);

  // and save the secondary momentum
  dp.sec = dp.mmom - dp.mom;
  dp.sec_beamcoord = dp.mmom_beamcoord - dp.mom_beamcoord;

  // The weight is the importance weight times the branching-ratio weight 
  weight = kaon.weight * br / SMKaonBR(kaon.meson_pdg);

  // set the mixing
  dp.C1 = fe;
  dp.C2 = 0;
  dp.C3 = 0.;
  dp.mass = fM;

  dp.meson_pdg = kaon.meson_pdg;
  dp.secondary_pdg = 22;
  dp.generator = 2; // kDP

  // equivalent neutrino energy
  dp.equiv_enu = EnuLab(flux.fnecm, dp.mmom, dp.pos);


  /*  // Get the DP Polarization 
  double meson_mass = Constants::Instance().kplus_mass;
  if(abs(dp.meson_pdg) == 321) {
    meson_mass = Constants::Instance().kplus_mass;
  }else if(abs(dp.meson_pdg) == 211) {
    meson_mass = Constants::Instance().piplus_mass;
  }
  */
  //  dp.polarization = PolDP(fM , lep_mass, meson_mass);
  //Account for the Polarization of m+ ->l+N being opposite to the polarization of m- -> l-bar{N} due to CP
  //  if(dp.secondary_pdg < 0) dp.polarization *= -1;

  return true;
}

DEFINE_ART_CLASS_TOOL(PiZero2DPFlux)

} // namespace ldm
} // namespace evgen
