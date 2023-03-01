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
#include "sbncode/EventGenerator/MeVPrtl/Tools/IMeVPrtlFlux.h"

#include "sbnobj/Common/EventGen/MeVPrtl/MeVPrtlFlux.h"


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
class MonoEnergyHiggsFlux : public IMeVPrtlFlux
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

    bool MakeFlux(const simb::MCFlux &flux, MeVPrtlFlux &higgs, double &weight) override;
    void configure(const fhicl::ParameterSet&) override;

    // no weights
    double MaxWeight() override { return -1.; }

private:
  TVector3 fStart; //!< Start of Higgs ray in detector coordinates [cm]
  TVector3 fDir; //!< Direction of Higgs ray (unit vector)
  double fE; //!< Energy of Higgs [GeV]
  double fM; //!< Mass of Higgs [GeV]
  double fMixingAngle;
  double fStartTime; //!< Start time of Higgs in detector time [us]
};

MonoEnergyHiggsFlux::MonoEnergyHiggsFlux(fhicl::ParameterSet const &pset):
  IMeVPrtlStage("MonoEnergyHiggsFlux"), 
  IMeVPrtlFlux(pset) 
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
  fStart = TVector3(pset.get<double>("X"), pset.get<double>("Y"), pset.get<double>("Z")); 
  fDir = TVector3(pset.get<double>("Xdir"), pset.get<double>("Ydir"), pset.get<double>("Zdir"));

  fE = pset.get<double>("E");
  fM = pset.get<double>("M");
  fMixingAngle = pset.get<double>("MixingAngle");
  fStartTime = pset.get<double>("T", 0.);

}

bool MonoEnergyHiggsFlux::MakeFlux(const simb::MCFlux &flux/*ignored*/, MeVPrtlFlux &higgs, double &weight) {
  (void) flux;
  higgs.pos = TLorentzVector(fStart, fStartTime);
  // set the momentum
  TVector3 p = fDir * sqrt(fE*fE - fM * fM);
  higgs.mom = TLorentzVector(p, fE);
  higgs.C1 = fMixingAngle;
  higgs.mass = fM;

  // no kaon here
  higgs.mmom = TLorentzVector(0, 0, 0, 0);
  higgs.meson_pdg = -1;

  // beam is same as detector coord
  higgs.pos_beamcoord = higgs.pos;
  higgs.mom_beamcoord = higgs.mom;
  higgs.mmom_beamcoord = higgs.mmom;
  higgs.generator = 0; // kDissonantHiggs


  weight = 1.;
  return true;
}

DEFINE_ART_CLASS_TOOL(MonoEnergyHiggsFlux)

} // namespace ldm
} // namespace evgen
