/**
 *
 */

// Framework Includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Utilities/ToolMacros.h"
#include "cetlib/cpu_timer.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// local includes
#include "IRayTrace.h"
#include "../Products/HiggsFlux.h"

// LArSoft includes
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/CoreUtils/ServiceUtil.h"

// std includes
#include <string>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace evgen {
namespace ldm {
/**
 *  @brief  RayTraceBox class definiton
 */
class RayTraceBox : public IRayTrace
{
public:
    /**
     *  @brief  Constructor
     */
    RayTraceBox(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~RayTraceBox();

    void configure(const fhicl::ParameterSet&) override;

    bool IntersectDetector(const HiggsFlux &flux, std::vector<TVector3> &intersection, float &weight) override;

    // no weights
    float ConstantWeight() override { return 1.; }
    float MaxWeight() override { return 1.; }

private:
  geo::BoxBoundedGeo fBox;
};

RayTraceBox::RayTraceBox(fhicl::ParameterSet const &pset):
  IHiggsStage("RayTraceBox") 
{
    this->configure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

RayTraceBox::~RayTraceBox()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
void RayTraceBox::configure(fhicl::ParameterSet const &pset)
{
  const geo::GeometryCore *geometry = lar::providerFrom<geo::Geometry>();
  fBox = geometry->DetectorEnclosureBox(pset.get<std::string>("Volume"));
}

    
bool RayTraceBox::IntersectDetector(const HiggsFlux &flux, std::vector<TVector3> &intersection, float &weight) {
  intersection = fBox.GetIntersections(flux.pos.Vect(), flux.mom.Vect().Unit());
  weight = 1.;
  return intersection.size() == 2;
}

DEFINE_ART_CLASS_TOOL(RayTraceBox)

} // namespace ldm
} // namespace evgen
