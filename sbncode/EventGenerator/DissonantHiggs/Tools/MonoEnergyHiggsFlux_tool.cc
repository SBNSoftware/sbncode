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

// local includes
#include "IHiggsFlux.h"
#include "../Products/HiggsFlux.h"

// LArSoft includes

// ROOT
#include "TVector3.h"

// std includes
#include <string>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace evgen {
namespace ldm {
/**
 *  @brief  MonoEnergyHiggsFlux class definiton
 */
class MonoEnergyHiggsFlux : public IHiggsFlux
{
public:
    /**
     *  @brief  Constructor
     */
    MonoEnergyHiggsFlux(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~MonoEnergyHiggsFlux();

    bool MakeFlux(const simb::MCFlux &flux, HiggsFlux &higgs, double &weight) override;
    void configure(const fhicl::ParameterSet&) override;

    // no weights
    float MaxWeight() override { return -1.; }

private:
  TVector3 fStart; //!< Start of Higgs ray in detector coordinates [cm]
  TVector3 fDir; //!< Direction of Higgs ray (unit vector)
  float fE; //!< Energy of Higgs [GeV]
  float fM; //!< Mass of Higgs [GeV]
  double fMixingAngle;
  float fStartTime; //!< Start time of Higgs in detector time [us]
};

MonoEnergyHiggsFlux::MonoEnergyHiggsFlux(fhicl::ParameterSet const &pset):
  IHiggsStage("MonoEnergyHiggsFlux") 
{
    this->configure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

MonoEnergyHiggsFlux::~MonoEnergyHiggsFlux()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
void MonoEnergyHiggsFlux::configure(fhicl::ParameterSet const &pset)
{
  fStart = TVector3(pset.get<float>("X"), pset.get<float>("Y"), pset.get<float>("Z")); 
  fDir = TVector3(pset.get<float>("Xdir"), pset.get<float>("Ydir"), pset.get<float>("Zdir"));

  fE = pset.get<float>("E");
  fM = pset.get<float>("M");
  fMixingAngle = pset.get<float>("MixingAngle");
  fStartTime = pset.get<float>("T", 0.);

}

bool MonoEnergyHiggsFlux::MakeFlux(const simb::MCFlux &flux/*ignored*/, HiggsFlux &higgs, double &weight) {
  (void) flux;
  higgs.pos = TLorentzVector(fStart, fStartTime);
  // set the momentum
  TVector3 p = fDir * sqrt(fE*fE - fM * fM);
  higgs.mom = TLorentzVector(p, fE);
  higgs.mixing = fMixingAngle;
  higgs.mass = fM;

  weight = 1.;
  return true;
}

DEFINE_ART_CLASS_TOOL(MonoEnergyHiggsFlux)

} // namespace ldm
} // namespace evgen
