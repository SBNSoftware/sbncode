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
#include "sbnobj/Common/EventGen/MeVPrtl/KaonParent.h"

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
 *
 *  Implementation of Kaon->HNL branching ratio taken from:
 *      arXiv:1912.07622
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

    bool MakeFlux(const simb::MCFlux &flux, MeVPrtlFlux &hnl, double &weight) override;
    void configure(const fhicl::ParameterSet&) override;

    double MaxWeight() override;

private:
  // config
  double fM; //!< Mass of HNL [GeV]
  double fMagUe4;
  double fMagUm4;
  double fTarget2Absorber;
  bool fKDAROnly;
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

// helper functions
//
// kaon -> leption + HNL
double hnl_momentum(double kaon_mass, double lep_mass, double hnl_mass) {
  if (kaon_mass - lep_mass < hnl_mass) return -1.;

  return sqrt(kaon_mass * kaon_mass * kaon_mass * kaon_mass 
    -2 * kaon_mass * kaon_mass * lep_mass * lep_mass
    -2 * kaon_mass * kaon_mass * hnl_mass * hnl_mass
       + lep_mass * lep_mass * lep_mass * lep_mass 
       + hnl_mass * hnl_mass * hnl_mass * hnl_mass 
    -2 * lep_mass * lep_mass * hnl_mass * hnl_mass) / ( 2 * kaon_mass );
}

double SMKaonBR(int kaon_pdg) {
  // The Kaons in Dk2nu file only include those that decay to neutrinos.
  //
  // We want all kaons -- in order to account for this, we divide by the 
  // branching-ratio of kaons to neutrinos
  //
  // Taken from: 
  // /cvmfs/minerva.opensciencegrid.org/minerva/beamsim/x86_64/geant4/source/particles/hadrons/mesons/src/G4KaonPlus.cc
  // /cvmfs/minerva.opensciencegrid.org/minerva/beamsim/x86_64/geant4/source/particles/hadrons/mesons/src/G4KaonZerLong.cc
  switch (kaon_pdg) {
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

double BranchingRatio(double hnl_mass, double u4, bool is_muon) {
  double lep_mass = is_muon ? Constants::Instance().muon_mass : Constants::Instance().elec_mass;
  double kplus_mass = Constants::Instance().kplus_mass;

  if (hnl_mass > kplus_mass - lep_mass) return 0.;

  double smbr = is_muon ? Constants::Instance().kaonp_mup_numu : Constants::Instance().kaonp_ep_nue;
  double lep_ratio = (lep_mass * lep_mass) / (kplus_mass * kplus_mass);
  double hnl_ratio = (hnl_mass * hnl_mass) / (kplus_mass * kplus_mass);
  double kinematic_factor = (lep_ratio + hnl_ratio - (lep_ratio - hnl_ratio) * (lep_ratio - hnl_ratio)) \
       * sqrt(1 + hnl_ratio * hnl_ratio +  lep_ratio * lep_ratio - 2*(hnl_ratio + lep_ratio + hnl_ratio*lep_ratio)) \
       / (lep_ratio * (1. - lep_ratio) * (1. - lep_ratio));

  // scale the branching ratio 
  return smbr * (u4 / (1. - u4)) * kinematic_factor;
}

std::pair<double, bool> Branch(double hnl_mass, double ue4, double um4, double rand) {
  double kplus_mass = Constants::Instance().kplus_mass;

  double br_muon = (um4 > 0. && hnl_mass < kplus_mass - Constants::Instance().muon_mass) ? BranchingRatio(hnl_mass, um4, true) : 0.;
  double br_elec = (ue4 > 0. && hnl_mass < kplus_mass - Constants::Instance().elec_mass) ? BranchingRatio(hnl_mass, ue4, false): 0.;
  if (br_muon == 0. && br_elec == 0.) return {0., false};

  double br_weight = br_muon + br_elec;

  bool is_muon = rand < (br_muon / br_weight);

  return {br_weight, is_muon};
}


//------------------------------------------------------------------------------------------------------------------------------------------
void Kaon2HNLFlux::configure(fhicl::ParameterSet const &pset)
{
  fM = pset.get<double>("M");
  fMagUm4 = pset.get<double>("MagUm4");
  fMagUe4 = pset.get<double>("MagUe4");

  fTarget2Absorber = pset.get<double>("Target2Absorber", 5000);
  fKDAROnly = pset.get<bool>("KDAROnly", false);

  double max_mass = (fMagUe4 > 0.) ? (Constants::Instance().kplus_mass - Constants::Instance().elec_mass) : 
      (Constants::Instance().kplus_mass - Constants::Instance().muon_mass);

  if (fM > max_mass) {
    throw cet::exception("Kaon2HNLFlux Tool: BAD MASS. Configured mass (" + std::to_string(fM) +
         ") is larger than maximum allowed by enabled couplings (" + std::to_string(max_mass) +  ").");
  }

}

double Kaon2HNLFlux::MaxWeight() { 
  // Weight comes from the NuMi importance weight -- max is 100 (add in an epsilon)
  // Scale by the branching ratios here
  return 100.0001 * std::max(BranchingRatio(fM, fMagUe4, false) / SMKaonBR(321), BranchingRatio(fM, fMagUm4, true) / SMKaonBR(321));
}

bool Kaon2HNLFlux::MakeFlux(const simb::MCFlux &flux, evgen::ldm::MeVPrtlFlux &hnl, double &weight) {
  // make the kaon parent
  evgen::ldm::KaonParent kaon(flux);
  if (abs(kaon.kaon_pdg) != 321) return false; // Only take charged kaons

  // select on the kaon
  if (fKDAROnly && (kaon.mom.P() > 1e-3 || kaon.pos.Z() < fTarget2Absorber)) return false;
  if (fKDAROnly) std::cout << "FOUND KDAR\n";

  TLorentzVector Beam4 = BeamOrigin();

  std::cout << "Beam Origin (x,y,z, toff): " << Beam4.X() 
					<< " " << Beam4.Y()
					<< " " << Beam4.Z()
					<< " " << Beam4.T()
  					<< std::endl;

  // get position in detector frame
  hnl.pos_beamcoord = kaon.pos;

  std::cout << "hnl pos beam coord: " << hnl.pos_beamcoord.X()
			<< " " << hnl.pos_beamcoord.Y() 
			<< " " << hnl.pos_beamcoord.Z() 
			<< " " << hnl.pos_beamcoord.T() 
			<< std::endl;

  hnl.pos = kaon.pos;
  hnl.pos.Transform(fBeam2Det);
  
  std::cout << "hnl pos det coord: " << hnl.pos.X()
			<< " " << hnl.pos.Y() 
			<< " " << hnl.pos.Z() 
			<< " " << hnl.pos.T() 
			<< std::endl;
  hnl.pos += Beam4;
  
  std::cout << "hnl pos det coord: " << hnl.pos.X()
			<< " " << hnl.pos.Y() 
			<< " " << hnl.pos.Z() 
			<< " " << hnl.pos.T() 
			<< std::endl;

  // Branch the parent Kaon Decay
  double hnl_mass = fM;
  std::pair<double, bool> decay = Branch(hnl_mass, fMagUe4, fMagUm4, GetRandom());
  double br = decay.first;
  bool is_muon = decay.second;
  double lep_mass = is_muon ? Constants::Instance().muon_mass : Constants::Instance().elec_mass;

  std::cout << "Kaon 2 HNL BR: " << br << std::endl;

  // ignore if we can't make this hnl
  // Ignore if branching ratio is exactly 0.
  if (br == 0.) return false;

  // get the momentum direction in the kaon parent rest frame
  double p = hnl_momentum(Constants::Instance().kplus_mass, lep_mass, hnl_mass);
  double e = sqrt(p*p + hnl_mass * hnl_mass);
  // Two-body decays are isotropic
  hnl.mom = TLorentzVector(p*RandomUnitVector(), e);

  // boost to lab frame
  TLorentzVector mom = hnl.mom;
  mom.Boost(kaon.mom.BoostVector());

  hnl.mom_beamcoord = mom;
  // rotate to detector frame
  hnl.mom = mom;
  hnl.mom.Transform(fBeam2Det);

  hnl.kmom_beamcoord = kaon.mom;
  // also save the kaon momentum in the detector frame
  hnl.kmom = kaon.mom;
  hnl.kmom.Transform(fBeam2Det);

  // and save the secondary momentum
  hnl.sec = hnl.kmom - hnl.mom;
  hnl.sec_beamcoord = hnl.kmom_beamcoord - hnl.mom_beamcoord;

  // The weight is the importance weight times the branching-ratio weight 
  weight = kaon.weight * br / SMKaonBR(kaon.kaon_pdg);

  // set the mixing
  hnl.C1 = fMagUe4;
  hnl.C2 = fMagUm4;
  hnl.mass = fM;

  hnl.kaon_pdg = kaon.kaon_pdg;
  hnl.secondary_pdg = (is_muon ? 11 : 13) * (kaon.kaon_pdg > 0 ? 1 : -1);
  hnl.generator = 1; // kHNL

  // equivalent neutrino energy
  hnl.equiv_enu = EnuLab(flux.fnecm, hnl.kmom, hnl.pos);

  return true;
}

DEFINE_ART_CLASS_TOOL(Kaon2HNLFlux)

} // namespace ldm
} // namespace evgen
