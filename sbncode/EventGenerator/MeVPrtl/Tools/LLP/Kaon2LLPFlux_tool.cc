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
 *  @brief  Kaon2LLPFlux class definiton
 *
 *  Implementation of Kaon->LLP branching ratio taken from:
 *      arXiv:1912.07622
 */
class Kaon2LLPFlux : public IMeVPrtlFlux
{
public:
    /**
     *  @brief  Constructor
     */
    Kaon2LLPFlux(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~Kaon2LLPFlux();

    bool MakeFlux(const simb::MCFlux &flux, MeVPrtlFlux &llp, double &weight) override;
    void configure(const fhicl::ParameterSet&) override;

    double MaxWeight() override;

private:
  // config
  double fM; //!< Mass of LLP [GeV]
  double fTarget2Absorber;
};

Kaon2LLPFlux::Kaon2LLPFlux(fhicl::ParameterSet const &pset):
  IMeVPrtlStage("Kaon2LLPFlux"), 
  IMeVPrtlFlux(pset) 
{
    this->configure(pset);

}

//------------------------------------------------------------------------------------------------------------------------------------------

Kaon2LLPFlux::~Kaon2LLPFlux()
{
}

// helper functions
//
// kaon -> lepton + LLP
double llp_momentum(double kaon_mass, double lep_mass, double llp_mass) {
  if (kaon_mass - lep_mass < llp_mass) return -1.;

  return sqrt(kaon_mass * kaon_mass * kaon_mass * kaon_mass 
    -2 * kaon_mass * kaon_mass * lep_mass * lep_mass
    -2 * kaon_mass * kaon_mass * llp_mass * llp_mass
       + lep_mass * lep_mass * lep_mass * lep_mass 
       + llp_mass * llp_mass * llp_mass * llp_mass 
    -2 * lep_mass * lep_mass * llp_mass * llp_mass) / ( 2 * kaon_mass );
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

//------------------------------------------------------------------------------------------------------------------------------------------
void Kaon2LLPFlux::configure(fhicl::ParameterSet const &pset)
{
  fM = pset.get<double>("M");

  fTarget2Absorber = pset.get<double>("Target2Absorber", 5000);

}

double Kaon2LLPFlux::MaxWeight() { 
  // Weight comes from the NuMi importance weight -- max is 100 (add in an epsilon)
  // Scale by the branching ratios here
  return 1;
}

bool Kaon2LLPFlux::MakeFlux(const simb::MCFlux &flux, evgen::ldm::MeVPrtlFlux &llp, double &weight) {
  // make the kaon parent
  evgen::ldm::MesonParent kaon(flux);
  if (abs(kaon.meson_pdg) != 321) return false; // Only take charged kaons

  TLorentzVector Beam4 = BeamOrigin();

  // get position in detector frame
  llp.pos_beamcoord = kaon.pos;
  llp.pos = kaon.pos;
  llp.pos.Transform(fBeam2Det);
  llp.pos += Beam4;

  // Branch the parent Kaon Decay
  double llp_mass = fM;
  double br = 1;
  bool is_muon = true;
  double lep_mass = is_muon ? Constants::Instance().muon_mass : Constants::Instance().elec_mass;

  if (fVerbose) std::cout << "BR: " << br << std::endl;

  // ignore if we can't make this llp
  // Ignore if branching ratio is exactly 0.
  if (br == 0.) return false;

  // get the momentum direction in the kaon parent rest frame
  double p = llp_momentum(Constants::Instance().kplus_mass, lep_mass, llp_mass);
  double e = sqrt(p*p + llp_mass * llp_mass);
  // Two-body decays are isotropic
  llp.mom = TLorentzVector(p*RandomUnitVector(), e);

  // boost to lab frame
  TLorentzVector mom = llp.mom;
  mom.Boost(kaon.mom.BoostVector());

  llp.mom_beamcoord = mom;
  // rotate to detector frame
  llp.mom = mom;
  llp.mom.Transform(fBeam2Det);

  llp.mmom_beamcoord = kaon.mom;
  // also save the kaon momentum in the detector frame
  llp.mmom = kaon.mom;
  llp.mmom.Transform(fBeam2Det);

  // and save the secondary momentum
  llp.sec = llp.mmom - llp.mom;
  llp.sec_beamcoord = llp.mmom_beamcoord - llp.mom_beamcoord;

  // The weight is the importance weight times the branching-ratio weight 
  weight = kaon.weight * br / SMKaonBR(kaon.meson_pdg);

  llp.mass = fM;

  llp.meson_pdg = kaon.meson_pdg;
  llp.secondary_pdg = (is_muon ? 13 : 11) * (kaon.meson_pdg > 0 ? 1 : -1);
  llp.generator = 1; // kLLP

  // equivalent neutrino energy
  llp.equiv_enu = EnuLab(flux.fnecm, llp.mmom, llp.pos);

  return true;

}
DEFINE_ART_CLASS_TOOL(Kaon2LLPFlux)

} // namespace ldm
} // namespace evgen
