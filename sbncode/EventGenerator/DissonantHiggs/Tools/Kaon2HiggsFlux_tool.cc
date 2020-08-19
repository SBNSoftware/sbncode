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

#define PION_MASS 139.57018 // MeV

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

    bool MakeFlux(const simb::MCFlux &flux, HiggsFlux &higgs, float &weight) override;
    void configure(const fhicl::ParameterSet&) override;

    // no weights
    float ConstantWeight() override { return 1.; }
    float MaxWeight() override { return 1.; }

private:
  // config
  float fM; //!< Mass of Higgs [GeV]
  std::string fSpillTimeConfig;

  // derived stuff
  evgb::EvtTimeShiftI *fTimeShiftMethod;
  TRotation fBeam2Det;
  TVector3 fBeamOrigin;
};

Kaon2HiggsFlux::Kaon2HiggsFlux(fhicl::ParameterSet const &pset):
  IHiggsStage("Kaon2HiggsFlux") 
{
    this->configure(pset);

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
}

//------------------------------------------------------------------------------------------------------------------------------------------

Kaon2HiggsFlux::~Kaon2HiggsFlux()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
void Kaon2HiggsFlux::configure(fhicl::ParameterSet const &pset)
{
  fM = pset.get<float>("M");
  fSpillTimeConfig = pset.get<std::string>("SpillTimeConfig", ""); 

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
}

bool Kaon2HiggsFlux::MakeFlux(const simb::MCFlux &flux, evgen::ldm::HiggsFlux &higgs, float &weight) {
  // make the kaon parent
  evgen::ldm::KaonParent kaon;
  bool success = evgen::ldm::MakeKaonParent(flux, kaon);
  if (!success) return false;

  float toff = fTimeShiftMethod ? fTimeShiftMethod->TimeOffset() : 0.;

  TLorentzVector Beam4(fBeamOrigin, toff);
  // get position in detector frame
  higgs.pos = kaon.pos.Transform(fBeam2Det) + Beam4;

  // get the momentum direction in the kaon parent rest frame
  //
  // isotropic decay
  float theta = CLHEP::RandFlat::shoot(fEngine, 0, M_PI);
  float phi = CLHEP::RandFlat::shoot(fEngine, 0, 2*M_PI);

  // calculate the momentum
  float kaon_mass = kaon.mom.M();  
  float higs_mass = fM;
  float pion_mass = PION_MASS;

  // ignore if we can't make this higgs
  if (kaon_mass - pion_mass < higs_mass) return false;

  float higgs_momentum = sqrt(kaon_mass * kaon_mass * kaon_mass * kaon_mass 
                         -2 * kaon_mass * kaon_mass * pion_mass * pion_mass
                         -2 * kaon_mass * kaon_mass * higs_mass * higs_mass
                            + pion_mass * pion_mass * pion_mass * pion_mass 
                            + higs_mass * higs_mass * higs_mass * higs_mass 
                         -2 * pion_mass * pion_mass * higs_mass * higs_mass) / ( 2 * kaon_mass );

  TVector3 kaon_frame_momentum;
  kaon_frame_momentum.SetMagThetaPhi(higgs_momentum, theta, phi);

  // boost to lab frame
  TLorentzVector mom;
  mom.SetVectM(kaon_frame_momentum, fM);
  mom.Boost(kaon.mom.BoostVector());

  // rotate to detector frame
  mom *= fBeam2Det;

  higgs.mom = mom;

  // no weighting here
  weight = 1;
  return true;
}

DEFINE_ART_CLASS_TOOL(Kaon2HiggsFlux)

#undef PION_MASS

} // namespace ldm
} // namespace evgen
