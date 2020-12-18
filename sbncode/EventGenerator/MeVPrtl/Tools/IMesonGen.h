/**
 *  @file   IMesonGen.h
 *
 *  @brief  This provides an art tool interface definition for tools which can create
 *          fake particles to overlay onto input daq fragments during decoding
 *
 *  @author grayputnam@uchicago.edu
 *
 */
#ifndef IMesonGen_h
#define IMesonGen_h

// Framework Includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"

#include "nusimdata/SimulationBase/MCFlux.h"

#include "IMeVPrtlStage.h"

// Algorithm includes

//------------------------------------------------------------------------------------------------------------------------------------------

namespace evgen
{
namespace ldm {
/**
 *  @brief  IMesonGen interface class definiton
 */
class IMesonGen: virtual public IMeVPrtlStage
{
public:
    /**
     *  @brief  Virtual Destructor
     */
    virtual ~IMesonGen() noexcept = default;

    virtual simb::MCFlux GetNext() = 0;
    virtual double GetPOT() = 0;

};

} // namespace ldm
} // namespace evgen
#endif

