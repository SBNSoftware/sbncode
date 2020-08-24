/**
 *  @file   IKaonGen.h
 *
 *  @brief  This provides an art tool interface definition for tools which can create
 *          fake particles to overlay onto input daq fragments during decoding
 *
 *  @author grayputnam@uchicago.edu
 *
 */
#ifndef IKaonGen_h
#define IKaonGen_h

// Framework Includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"

#include "nusimdata/SimulationBase/MCFlux.h"

#include "IHiggsStage.h"

// Algorithm includes

//------------------------------------------------------------------------------------------------------------------------------------------

namespace evgen
{
namespace ldm {
/**
 *  @brief  IKaonGen interface class definiton
 */
class IKaonGen: virtual public IHiggsStage
{
public:
    /**
     *  @brief  Virtual Destructor
     */
    virtual ~IKaonGen() noexcept = default;

    virtual simb::MCFlux GetNext() = 0;
    virtual double GetPOT() = 0;

};

} // namespace ldm
} // namespace evgen
#endif

