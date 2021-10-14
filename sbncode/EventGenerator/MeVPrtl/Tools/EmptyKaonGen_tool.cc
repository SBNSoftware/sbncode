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

#include "nugen/EventGeneratorBase/GENIE/EvtTimeShiftFactory.h"
#include "nugen/EventGeneratorBase/GENIE/EvtTimeShiftI.h"
#include "nusimdata/SimulationBase/MCFlux.h"

// local includes
#include "IMesonGen.h"
#include "sbnobj/Common/EventGen/MeVPrtl/MeVPrtlFlux.h"

// LArSoft includes
#include "dk2nu/tree/dk2nu.h"
#include "dk2nu/tree/dkmeta.h"
#include "dk2nu/tree/NuChoice.h"

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
 *  @brief  EmptyKaonGen class definiton
 */
class EmptyKaonGen : public IMesonGen
{
public:
    /**
     *  @brief  Constructor
     */
    EmptyKaonGen(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~EmptyKaonGen();

    simb::MCFlux GetNext() override;
    void configure(const fhicl::ParameterSet&) override;

    // no POT
    double GetPOT() override { return 0.; }

    // no weights
    double MaxWeight() override  { return -1.; }
};

EmptyKaonGen::EmptyKaonGen(fhicl::ParameterSet const &pset):
  IMeVPrtlStage("EmptyKaonGen") 
{}

//------------------------------------------------------------------------------------------------------------------------------------------

EmptyKaonGen::~EmptyKaonGen() {}

//------------------------------------------------------------------------------------------------------------------------------------------
void EmptyKaonGen::configure(fhicl::ParameterSet const &pset) {}

simb::MCFlux EmptyKaonGen::GetNext() {
  simb::MCFlux flux;
  return flux;
}

DEFINE_ART_CLASS_TOOL(EmptyKaonGen)

} // namespace ldm
} // namespace evgen
