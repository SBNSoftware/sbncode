/**
 *  @file   IMesonGen.h
 *
 *  @brief  This is an interface for an art Tool which sources MCFlux objects for
 *  downstream processing and tabulates POT information.
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

