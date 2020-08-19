/**
 *  @file   IHiggsDecay.h
 *
 *  @brief  This provides an art tool interface definition for tools which can create
 *          fake particles to overlay onto input daq fragments during decoding
 *
 *  @author grayputnam@uchicago.edu
 *
 */
#ifndef IHiggsDecay_h
#define IHiggsDecay_h

// Framework Includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"

#include "../Products/HiggsFlux.h"

#include "IHiggsStage.h"

// Algorithm includes
#include "nusimdata/SimulationBase/MCTruth.h"

#include <utility>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace evgen
{
namespace ldm {
/**
 *  @brief  IHiggsDecay interface class definiton
 */
class IHiggsDecay: virtual public IHiggsStage
{
public:
    /**
     *  @brief  Virtual Destructor
     */
    virtual ~IHiggsDecay() noexcept = default;

    virtual bool Decay(const HiggsFlux &flux, const TVector3 &in, const TVector3 &out, simb::MCTruth &truth, float &weight) = 0;
};

} // namespace ldm
} // namespace evgen
#endif

