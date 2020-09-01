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
#include "IHiggsFlux.h"
#include "../Products/HiggsFlux.h"
#include "../Products/KaonParent.h"

// LArSoft includes
#include "nugen/EventGeneratorBase/GENIE/EvtTimeShiftFactory.h"
#include "nugen/EventGeneratorBase/GENIE/EvtTimeShiftI.h"

// ROOT
#include "TVector3.h"

// std includes
#include <string>
#include <iostream>
#include <memory>

#include "../ParticleData.h"

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace evgen {
namespace ldm {
/**
 *  @brief  Kaon2HiggsFlux class definiton
 */
class Kaon2HiggsFlux : public IHiggsFlux
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

    bool MakeFlux(const simb::MCFlux &flux, HiggsFlux &higgs, double &weight) override;
    void configure(const fhicl::ParameterSet&) override;

    float MaxWeight() override { 
      // Weight comes from the NuMi importance weight -- max is 100 
      //
      // Also get the max BR
      return 100. * std::max(fKPBR, fKLBR);
    }

private:
  // config
  float fM; //!< Mass of Higgs [GeV]
  float fMixingAngle; //!< Mixing angle of dark higgs
  std::string fSpillTimeConfig;
  bool fKDAROnly;

  // derived stuff
  evgb::EvtTimeShiftI *fTimeShiftMethod;
  TRotation fBeam2Det;
  TVector3 fBeamOrigin;

  // branching ratios
  double fKPBR;
  double fKLBR;
};

Kaon2HiggsFlux::Kaon2HiggsFlux(fhicl::ParameterSet const &pset):
  IHiggsStage("Kaon2HiggsFlux") 
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

static const double tquark_mass = 172.76; // GeV
static const double higgs_vev = 246.22; // GeV

// CKM matrix
static const double abs_VtsVtd_squared = 1.0185e-07;
static const double rel_VtsVtd_squared = 1.0185e-07;

static const double hbar = 6.582119569e-25; // GeV*s
static const double c_cm_per_ns = 29.9792; // cm / ns

// compute branching ratios
double KaonPlusBranchingRatio(double higs_mass, double mixing) {
  // Kplus width
  //
  // matrix element for kplus
  double M_KP = (1. / 2.) * ( 3. / (16. * M_PI * M_PI * higgs_vev * higgs_vev * higgs_vev)) * (kplus_mass * kplus_mass) * (tquark_mass * tquark_mass) * mixing;
  double M_KP2 = M_KP * M_KP * abs_VtsVtd_squared;

  double kplus_width = (2 * higgs_momentum(kplus_mass, pplus_mass, higs_mass)/kplus_mass) * M_KP2 / (16 * M_PI * kplus_mass);

  // convert to partial lifetime
  double kplus_2higgs_lifetime = hbar / kplus_width; 

  // and branching ratio
  //
  // (this higgs decay would make a negligible contribution to the overall lifetime)
  return kplus_lifetime / kplus_2higgs_lifetime;
}

double KaonLongBranchingRatio(double higs_mass, double mixing) {
  double M_KL = (1. / 2.) * (3. / (16. * M_PI * M_PI * higgs_vev * higgs_vev * higgs_vev)) * (klong_mass * klong_mass) * (tquark_mass * tquark_mass) * mixing;
  double M_KL2 = M_KL * M_KL * rel_VtsVtd_squared;

  double klong_width = (2 * higgs_momentum(klong_mass, pzero_mass, higs_mass) / klong_mass) * M_KL2 / (16 * M_PI * klong_mass);

  double klong_2higgs_lifetime = hbar / klong_width;

  return klong_lifetime / klong_2higgs_lifetime;

}


//------------------------------------------------------------------------------------------------------------------------------------------
void Kaon2HiggsFlux::configure(fhicl::ParameterSet const &pset)
{
  fM = pset.get<float>("M");
  fMixingAngle = pset.get<float>("MixingAngle");
  fSpillTimeConfig = pset.get<std::string>("SpillTimeConfig", ""); 

  fKDAROnly = pset.get<bool>("KDAROnly", false);

  // rotation matrix
  std::vector<double> rotation = pset.get<std::vector<double>>("Beam2DetectorRotation");
  assert(rotation.size() == 9);
  // get the three axes
  TVector3 beam2DetX(rotation[0], rotation[3], rotation[6]);
  TVector3 beam2DetY(rotation[1], rotation[4], rotation[7]);
  TVector3 beam2DetZ(rotation[2], rotation[5], rotation[8]);

  fBeam2Det.RotateAxes(beam2DetX, beam2DetY, beam2DetZ);

  // beam origin
  std::vector<double> origin = pset.get<std::vector<double>>("BeamOrigin");
  assert(origin.size() == 3 || origin.size() == 6);
  // beam origin specified in detector rotation-frame
  if (origin.size() == 3) {
    fBeamOrigin.SetXYZ(origin[0], origin[1], origin[2]);
  }
  // beam origin specified in beam rotation-frame
  else if (origin.size() == 6) {
    TVector3 userpos = TVector3(origin[0], origin[1], origin[2]);
    TVector3 beampos = TVector3(origin[3], origin[4], origin[5]);
    fBeamOrigin = userpos - fBeam2Det * beampos;
  }

  fKPBR = KaonPlusBranchingRatio(fM, fMixingAngle);
  fKLBR = KaonLongBranchingRatio(fM, fMixingAngle);

  // get random seed for stuff
  unsigned seed = art::ServiceHandle<rndm::NuRandomService>()->getSeed();

  // use the time-shifting tools from GENIE
  fTimeShiftMethod = NULL;
  if (fSpillTimeConfig != "") {
    fTimeShiftMethod = evgb::EvtTimeShiftFactory::Instance().GetEvtTimeShift(fSpillTimeConfig);
    if ( fTimeShiftMethod ) {
      if ( ! fTimeShiftMethod->IsRandomGeneratorSeeded() ) {
        fTimeShiftMethod->GetRandomGenerator()->SetSeed(seed);
      }
      fTimeShiftMethod->PrintConfig();
    } 
    else {
      evgb::EvtTimeShiftFactory::Instance().Print();
    }
  }

  std::cout << "K+ branching ratio: " << fKPBR << std::endl;
  std::cout << "K0 branching ratio: " << fKLBR << std::endl;

  std::cout << "Neutrino TIF: " << (fBeamOrigin.Mag()/c_cm_per_ns) << std::endl;

  if (fTimeShiftMethod) {
    std::cout << "Timing Config:\n";
    fTimeShiftMethod->PrintConfig();
    std::cout << std::endl;
  }
}

bool Kaon2HiggsFlux::MakeFlux(const simb::MCFlux &flux, evgen::ldm::HiggsFlux &higgs, double &weight) {
  // make the kaon parent
  evgen::ldm::KaonParent kaon;
  bool success = evgen::ldm::MakeKaonParent(flux, kaon);
  if (!success) return false;

  // select on the kaon
  if (fKDAROnly && (kaon.mom.P() > 1e-3 || kaon.pos.Z() < 72000.)) return false;
  if (fKDAROnly) std::cout << "FOUND KDAR\n";

  float toff = fTimeShiftMethod ? fTimeShiftMethod->TimeOffset() : 0.;

  // subtract out the delay of neutrinos reaching the beam
  float neutrino_tif = fBeamOrigin.Mag()/c_cm_per_ns;
  toff -= neutrino_tif;

  TLorentzVector Beam4(fBeamOrigin, toff);
  // get position in detector frame
  higgs.pos_beamcoord = kaon.pos;
  higgs.pos = kaon.pos;
  higgs.pos.Transform(fBeam2Det);
  higgs.pos += Beam4;

  // get the momentum direction in the kaon parent rest frame
  float kaon_mass = kaon.mom.M();  
  float higs_mass = fM;
  float pion_mass = evgen::ldm::PDATA->GetParticle(kaon.pion_pdg)->Mass();

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

  higgs.kmom_beamcoord = kaon.mom;
  // also save the kaon momentum in the detector frame
  higgs.kmom = kaon.mom;
  higgs.kmom.Transform(fBeam2Det);

  // The weight is the importance weight times the branching-ratio weight 
  weight = kaon.weight;
  if (kaon.kaon_pdg == 130 /* KLong */) {
    weight = weight * fKLBR;
  }
  else { // KPlus or KMinus
    weight = weight * fKPBR;
  }

  // set the mixing
  higgs.mixing = fMixingAngle;
  higgs.mass = fM;

  higgs.kaon_pdg = kaon.kaon_pdg;

  return true;
}

DEFINE_ART_CLASS_TOOL(Kaon2HiggsFlux)

} // namespace ldm
} // namespace evgen
