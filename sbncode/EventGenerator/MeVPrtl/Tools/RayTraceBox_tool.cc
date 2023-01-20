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
#include "sbnobj/Common/EventGen/MeVPrtl/MeVPrtlFlux.h"

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

    bool IntersectDetector(MeVPrtlFlux &flux, std::array<TVector3, 2> &intersection, double &weight) override;

    // no weights
    double MaxWeight() override { return 1.; }

private:
  geo::BoxBoundedGeo fBox;
  bool fVerbose;
};

RayTraceBox::RayTraceBox(fhicl::ParameterSet const &pset):
  IMeVPrtlStage("RayTraceBox") 
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
  fVerbose = pset.get<bool>("Verbose", true);  
  
  if (pset.has_key("Box")) {
    std::array<double, 6> box_config = pset.get<std::array<double, 6>>("Box");
    // xmin, xmax, ymin, ymax, zmin, zmax
    fBox = geo::BoxBoundedGeo(box_config[0], box_config[1], box_config[2], box_config[3], box_config[4], box_config[5]);
  }
  else {
    const geo::GeometryCore *geometry = lar::providerFrom<geo::Geometry>();
    fBox = geometry->DetectorEnclosureBox(pset.get<std::string>("Volume"));
  }

  if (fVerbose){
    std::cout << "Detector Box." << std::endl;
    std::cout << "X " << fBox.MinX() << " " << fBox.MaxX() << std::endl;
    std::cout << "Y " << fBox.MinY() << " " << fBox.MaxY() << std::endl;
    std::cout << "Z " << fBox.MinZ() << " " << fBox.MaxZ() << std::endl;
  }
}

    
bool RayTraceBox::IntersectDetector(MeVPrtlFlux &flux, std::array<TVector3, 2> &intersection, double &weight) {
  std::vector<TVector3> box_intersections = fBox.GetIntersections(flux.pos.Vect(), flux.mom.Vect().Unit());

  if (box_intersections.size() != 2) return false;

  TVector3 A = box_intersections[0];
  TVector3 B = box_intersections[1];

  // make sure that the flux start lies outside the detector
  if ((flux.pos.Vect() - A).Mag() < (A-B).Mag() && (flux.pos.Vect() - B).Mag() < (A-B).Mag()) {
    throw cet::exception("RayTraceBox Exception", "Input portal flux starts inside detector volume: "
        "MeVPrtl start At: (" + std::to_string(flux.pos.X()) + ", " + std::to_string(flux.pos.Y()) + ", " + std::to_string(flux.pos.Z()) + "). "
        "Intersection A At: " + std::to_string(A.X()) + ", " + std::to_string(A.Y()) + ", " + std::to_string(A.Z()) + "). "
        "Intersection B At: " + std::to_string(B.X()) + ", " + std::to_string(B.Y()) + ", " + std::to_string(B.Z()) + ").\n"
    );
  } 

  // if the ray points the wrong way, it doesn't intersect
  if (flux.mom.Vect().Unit().Dot((A - flux.pos.Vect()).Unit()) < 0.) {
    if (fVerbose) std::cout << "RAYTRACE: MeVPrtl points wrong way" << std::endl;
    return false;
  }

  if ((flux.pos.Vect() - A).Mag() < (flux.pos.Vect() - B).Mag()) {
    intersection = {A, B}; // A is entry, B is exit
  }
  else {
    intersection = {B, A}; // reversed
  }

  weight = 1.;
  return true;
}

DEFINE_ART_CLASS_TOOL(RayTraceBox)

} // namespace ldm
} // namespace evgen
