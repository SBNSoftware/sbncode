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
 *  @brief  Kaon2HiggsFlux class definiton
 *  Implementation of model taken from:
 *      https://arxiv.org/abs/1909.11670
 */
class Kaon2HiggsFlux : public IMeVPrtlFlux
{
public:
    /**
     *  @brief  Constructor
     */
    Kaon2HiggsFlux(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~Kaon2HiggsFlux();

    bool MakeFlux(const simb::MCFlux &flux, MeVPrtlFlux &higgs, double &weight) override;
    void configure(const fhicl::ParameterSet&) override;

    double MaxWeight() override; 

private:
  // config
  double fM; //!< Mass of Higgs [GeV]
  double fMixingAngle; //!< Mixing angle of dark higgs
  bool fKDAROnly;
  bool fKDIFOnly;
  bool fKDIFandBeamline;
  bool fIgnoreParentDecayTime;

  // branching ratios
  double fKPBR;
  double fKLBR;
};

Kaon2HiggsFlux::Kaon2HiggsFlux(fhicl::ParameterSet const &pset):
  IMeVPrtlStage("Kaon2HiggsFlux"), 
  IMeVPrtlFlux(pset) 
{
    this->configure(pset);

}

//------------------------------------------------------------------------------------------------------------------------------------------

Kaon2HiggsFlux::~Kaon2HiggsFlux()
{
}

double higgs_momentum(double kaon_mass, double pion_mass, double higs_mass) {
  return sqrt(kaon_mass * kaon_mass * kaon_mass * kaon_mass 
    -2 * kaon_mass * kaon_mass * pion_mass * pion_mass
    -2 * kaon_mass * kaon_mass * higs_mass * higs_mass
       + pion_mass * pion_mass * pion_mass * pion_mass 
       + higs_mass * higs_mass * higs_mass * higs_mass 
    -2 * pion_mass * pion_mass * higs_mass * higs_mass) / ( 2 * kaon_mass );
}

// compute branching ratios
double KaonPlusBranchingRatio(double higs_mass, double mixing) {
  double kplus_mass = Constants::Instance().kplus_mass;
  double piplus_mass = Constants::Instance().piplus_mass;
  double tquark_mass = Constants::Instance().tquark_mass;
  double higgs_vev = Constants::Instance().higgs_vev;
  double abs_VtsVtd_squared = Constants::Instance().abs_VtsVtd_squared;

  // Kplus width
  //
  // matrix element for kplus
  double M_KP = (1. / 2.) * ( 3. / (16. * M_PI * M_PI * higgs_vev * higgs_vev * higgs_vev)) * (kplus_mass * kplus_mass) * (tquark_mass * tquark_mass) * mixing;
  double M_KP2 = M_KP * M_KP * abs_VtsVtd_squared;

  double kplus_width = (2 * higgs_momentum(kplus_mass, piplus_mass, higs_mass)/kplus_mass) * M_KP2 / (16 * M_PI * kplus_mass);

  // convert to partial lifetime
  double kplus_2higgs_lifetime = Constants::Instance().hbar / kplus_width; 

  // and branching ratio
  //
  // (this higgs decay would make a negligible contribution to the overall lifetime)
  return Constants::Instance().kplus_lifetime / kplus_2higgs_lifetime;
}

double KaonLongBranchingRatio(double higs_mass, double mixing) {
  double klong_mass = Constants::Instance().klong_mass;
  double pizero_mass = Constants::Instance().pizero_mass;
  double tquark_mass = Constants::Instance().tquark_mass;
  double higgs_vev = Constants::Instance().higgs_vev;
  double rel_VtsVtd_squared = Constants::Instance().rel_VtsVtd_squared;

  double M_KL = (1. / 2.) * (3. / (16. * M_PI * M_PI * higgs_vev * higgs_vev * higgs_vev)) * (klong_mass * klong_mass) * (tquark_mass * tquark_mass) * mixing;
  double M_KL2 = M_KL * M_KL * rel_VtsVtd_squared;

  double klong_width = (2 * higgs_momentum(klong_mass, pizero_mass, higs_mass) / klong_mass) * M_KL2 / (16 * M_PI * klong_mass);

  double klong_2higgs_lifetime = Constants::Instance().hbar / klong_width;

  return Constants::Instance().klong_lifetime / klong_2higgs_lifetime;

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

int PionPdg(int kaon_pdg) {
  switch (kaon_pdg) {
    case 321:
      return 211;
    case -321:
      return -211;
    case 130:
      return 111;
    default: 
      return -1;
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------
void Kaon2HiggsFlux::configure(fhicl::ParameterSet const &pset)
{
  fM = pset.get<double>("M");
  fMixingAngle = pset.get<double>("MixingAngle");
  fIgnoreParentDecayTime = pset.get<bool>("IgnoreParentDecayTime");
  fKDAROnly = pset.get<bool>("KDAROnly", false);
  fKDIFOnly = pset.get<bool>("KDIFOnly", false);
  fKDIFandBeamline = pset.get<bool>("KDIFandBeamline", false);

  // Throw exception for a bad mass value
  if (fM > Constants::Instance().kplus_mass - Constants::Instance().piplus_mass && 
      fM > Constants::Instance().klong_mass - Constants::Instance().pizero_mass) {
    throw cet::exception("Kaon2HiggsFlux Tool: BAD MASS. Configured mass (" + std::to_string(fM) +
         ") is larger than allowed for K+ (" + std::to_string(Constants::Instance().kplus_mass - Constants::Instance().piplus_mass) + ") and KL (" +
         std::to_string(Constants::Instance().klong_mass - Constants::Instance().pizero_mass) + " production.");
  } 

  fKPBR = KaonPlusBranchingRatio(fM, fMixingAngle);
  fKLBR = KaonLongBranchingRatio(fM, fMixingAngle);

  std::cout << "K+ branching ratio: " << fKPBR << std::endl;
  std::cout << "K0 branching ratio: " << fKLBR << std::endl;

}

double Kaon2HiggsFlux::MaxWeight() {
  // Weight comes from the NuMi importance weight -- max is 100 (add in an epsilon) 
  //
  // Also get the max BR
  return 100.0001 * std::max(fKPBR / SMKaonBR(321), fKLBR / SMKaonBR(130));
}


bool Kaon2HiggsFlux::MakeFlux(const simb::MCFlux &flux, evgen::ldm::MeVPrtlFlux &higgs, double &weight) {

  // make the kaon parent
  evgen::ldm::MesonParent kaon(flux);
  if (!kaon.isKaon()) return false; // parent wasn't a kaon

  // select on the kaon
  if (fKDAROnly && (kaon.mom.P() > 1e-3 || kaon.pos.Z() < 72000.)) return false; //selects KDAR from absorber only.
  if (fKDAROnly) std::cout << "FOUND KDAR\n";
  if (fKDIFOnly && (kaon.mom.P() <= 1e-3)){ // no KDAR allowed (from anywhere). Accepts KDIF from beamline or absorber.
    std::cout << "found KDAR, skipping to next event\n";
    return false;
  }
  if (fKDIFandBeamline && (kaon.mom.P() <= 1e-3 && kaon.pos.Z() >= 72000.)) return false; //allows for KDAR from decay-pipe, and any KDIF. This option is exactly orthogonal to "KDAROnly" option. 

  TLorentzVector Beam4 = BeamOrigin();
  // get position in detector frame
  higgs.pos_beamcoord = kaon.pos;
  higgs.pos = kaon.pos;
  higgs.pos.Transform(fBeam2Det);
  higgs.pos += Beam4;

  if (fIgnoreParentDecayTime) higgs.pos.SetT(Beam4.T());

  // get the momentum direction in the kaon parent rest frame
  double kaon_mass = kaon.mom.M();  
  double higs_mass = fM;
  double pion_mass = TDatabasePDG::Instance()->GetParticle(PionPdg(kaon.meson_pdg))->Mass();

  // ignore if we can't make this higgs
  if (kaon_mass - pion_mass < higs_mass) return false;

  // isotropic decay
  TVector3 kaon_frame_momentum = RandomUnitVector() * higgs_momentum(kaon_mass, pion_mass, higs_mass);
  std::cout << "Rest frame higgs P: " <<  higgs_momentum(kaon_mass, pion_mass, higs_mass) << std::endl;

  // boost to lab frame
  TLorentzVector mom;
  mom.SetVectM(kaon_frame_momentum, higs_mass);
  mom.Boost(kaon.mom.BoostVector());

  higgs.mom_beamcoord = mom;
  // rotate to detector frame
  higgs.mom = mom;
  higgs.mom.Transform(fBeam2Det);

  higgs.mmom_beamcoord = kaon.mom;
  // also save the kaon momentum in the detector frame
  higgs.mmom = kaon.mom;
  higgs.mmom.Transform(fBeam2Det);

  // The weight is the importance weight times the branching-ratio weight 
  weight = kaon.weight / SMKaonBR(kaon.meson_pdg);
  if (kaon.meson_pdg == 130 /* KLong */) {
    weight = weight * fKLBR;
  }
  else { // KPlus or KMinus
    weight = weight * fKPBR;
  }

  // and save the secondary momentum
  higgs.sec = higgs.mmom - higgs.mom;
  higgs.sec_beamcoord = higgs.mmom_beamcoord - higgs.mom_beamcoord;

  // set the mixing
  higgs.C1 = fMixingAngle;
  higgs.mass = fM;

  higgs.meson_pdg = kaon.meson_pdg;
  higgs.generator = 0; // kDissonantHiggs
  higgs.secondary_pdg = PionPdg(kaon.meson_pdg);
  higgs.equiv_enu = EnuLab(flux.fnecm, higgs.mmom, higgs.pos);

  return true;
}

DEFINE_ART_CLASS_TOOL(Kaon2HiggsFlux)

} // namespace ldm
} // namespace evgen
