/**
 *  @file   IHiggsFlux.h
 *
 *  @brief  This provides an art tool interface definition for tools which can create
 *          fake particles to overlay onto input daq fragments during decoding
 *
 *  @author grayputnam@uchicago.edu
 *
 */
#ifndef IHiggsFlux_h
#define IHiggsFlux_h

// Framework Includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"

#include "../Products/HiggsFlux.h"
#include "nusimdata/SimulationBase/MCFlux.h"

#include "IHiggsStage.h"

// Algorithm includes

//------------------------------------------------------------------------------------------------------------------------------------------

namespace evgen
{
namespace ldm {
/**
 *  @brief  IHiggsFlux interface class definiton
 */
class IHiggsFlux: virtual public IHiggsStage
{
public:
    /**
     *  @brief  Virtual Destructor
     */
    virtual ~IHiggsFlux() noexcept = default;

    virtual bool MakeFlux(const simb::MCFlux &flux, HiggsFlux &higgs, float &weight) = 0;

};

} // namespace ldm
} // namespace evgen
#endif

