/**
 *  @file   IRayTrace.h
 *
 *  @brief  This provides an interface for an art tool which ray traces "Prtl"
 *  (massive) particles from their production position to impinge on a detector
 *  surface. It provides the intersection points of a MeVPrtlFlux object, along with
 *  a weight. 
 *
 *  @author grayputnam@uchicago.edu
 *
 */
#ifndef IRayTrace_h
#define IRayTrace_h

// Framework Includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"

// Algorithm includes
#include "TVector3.h"

#include "IMeVPrtlStage.h"

// cpp includes
#include <vector>
#include <array>

#include "../Products/MeVPrtlFlux.h"

//------------------------------------------------------------------------------------------------------------------------------------------

namespace evgen
{
namespace ldm {
/**
 *  @brief  IRayTrace interface class definiton
 */
class IRayTrace: virtual public IMeVPrtlStage
{
public:
    /**
     *  @brief  Virtual Destructor
     */
    virtual ~IRayTrace() noexcept = default;

    virtual bool IntersectDetector(MeVPrtlFlux &flux, std::array<TVector3, 2> &intersect, double &weight) = 0;
};

} // namespace ldm
} // namespace evgen
#endif

