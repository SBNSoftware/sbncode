/**
 *  @file   IMeVPrtlDecay.h
 *
 *  @brief  This provides an art tool interface definition for tools which can create
 *          fake particles to overlay onto input daq fragments during decoding
 *
 *  @author grayputnam@uchicago.edu
 *
 */
#ifndef IMeVPrtlDecay_h
#define IMeVPrtlDecay_h

// Framework Includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"

#include "../Products/MeVPrtlFlux.h"
#include "../Products/MeVPrtlDecay.h"

#include "IMeVPrtlStage.h"
#include "Constants.h"

// Algorithm includes

#include <utility>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace evgen
{
namespace ldm {

/**
 *  @brief  IMeVPrtlDecay interface class definiton
 */
class IMeVPrtlDecay: virtual public IMeVPrtlStage
{
public:
    /**
     *  @brief  Virtual Destructor
     */
    virtual ~IMeVPrtlDecay() noexcept = default;

    virtual bool Decay(const MeVPrtlFlux &flux, const TVector3 &in, const TVector3 &out, MeVPrtlDecay &decay, double &weight) = 0;

protected:
    double TimeOfFlight(const MeVPrtlFlux &flux, TVector3 decay) {
      // TODO: should the neutrino TOF be subtracted here to get the correct T0?
      return flux.pos.T() + (flux.pos.Vect() - decay).Mag() * (1. / flux.mom.Beta()) / c_cm_per_ns;
    }
};

} // namespace ldm
} // namespace evgen
#endif

