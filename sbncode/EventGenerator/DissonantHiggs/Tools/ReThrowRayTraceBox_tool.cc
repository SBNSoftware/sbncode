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
 *  @brief  ReThrowRayTraceBox class definiton
 */
class ReThrowRayTraceBox : public IRayTrace
{
public:
    /**
     *  @brief  Constructor
     */
    ReThrowRayTraceBox(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~ReThrowRayTraceBox();

    void configure(const fhicl::ParameterSet&) override;

    bool IntersectDetector(HiggsFlux &flux, std::array<TVector3, 2> &intersection, double &weight) override;

    TLorentzVector ThrowHiggsMomentum(const HiggsFlux &flux);

    // always thrown at least once
    float MaxWeight() override { return 1.; }

private:
  geo::BoxBoundedGeo fBox;
  unsigned fNThrows;
};

ReThrowRayTraceBox::ReThrowRayTraceBox(fhicl::ParameterSet const &pset):
  IHiggsStage("ReThrowRayTraceBox") 
{
    this->configure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

ReThrowRayTraceBox::~ReThrowRayTraceBox()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
void ReThrowRayTraceBox::configure(fhicl::ParameterSet const &pset)
{
  const geo::GeometryCore *geometry = lar::providerFrom<geo::Geometry>();
  fBox = geometry->DetectorEnclosureBox(pset.get<std::string>("Volume"));

  fNThrows = pset.get<unsigned>("NThrows", 100000);

  std::cout << "Detector Box." << std::endl;
  std::cout << "X " << fBox.MinX() << " " << fBox.MaxX() << std::endl;
  std::cout << "Y " << fBox.MinY() << " " << fBox.MaxY() << std::endl;
  std::cout << "Z " << fBox.MinZ() << " " << fBox.MaxZ() << std::endl;

}

TLorentzVector ReThrowRayTraceBox::ThrowHiggsMomentum(const HiggsFlux &flux) {
  // pick a direction for the rest-frame
  TVector3 dir = RandomUnitVector();

  // make the higgs momentum this
  // 
  // Move the higgs momentum back to the kaon rest frame
  TLorentzVector higgs_mom = flux.mom;
  higgs_mom.Boost(-flux.kmom.BoostVector());
  
  higgs_mom = TLorentzVector(higgs_mom.P() * dir, higgs_mom.E());

  // boost back
  higgs_mom.Boost(flux.kmom.BoostVector());

  return higgs_mom;
  
}
    
bool ReThrowRayTraceBox::IntersectDetector(HiggsFlux &flux, std::array<TVector3, 2> &intersection, double &weight) {
  // try out the higgs direction a bunch of times
  unsigned i = 1;
  std::vector<TVector3> box_intersections;
  TLorentzVector higgs_mom;
  for (; i <= fNThrows; i++) {
    higgs_mom = ThrowHiggsMomentum(flux);
    box_intersections = fBox.GetIntersections(flux.pos.Vect(), higgs_mom.Vect().Unit());
    if (box_intersections.size() == 2) {
      break;
    }
  }

  // did we get a hit?
  if (box_intersections.size() != 2) return false;

  TVector3 A = box_intersections[0];
  TVector3 B = box_intersections[1];

  // make sure that the flux start lies outside the detector
  if ((flux.pos.Vect() - A).Mag() < (A-B).Mag() && (flux.pos.Vect() - B).Mag() < (A-B).Mag()) {
    throw cet::exception("ReThrowRayTraceBox Exception", "Input higgs flux starts inside detector volume: "
        "Higgs start At: (" + std::to_string(flux.pos.X()) + ", " + std::to_string(flux.pos.Y()) + ", " + std::to_string(flux.pos.Z()) + "). "
        "Intersection A At: " + std::to_string(A.X()) + ", " + std::to_string(A.Y()) + ", " + std::to_string(A.Z()) + "). "
        "Intersection B At: " + std::to_string(B.X()) + ", " + std::to_string(B.Y()) + ", " + std::to_string(B.Z()) + ").\n"
    );
  } 

  // if the ray points the wrong way, it doesn't intersect
  if (flux.mom.Vect().Unit().Dot((A - flux.pos.Vect()).Unit()) < 0.) {
    std::cout << "RAYTRACE: Higgs points wrong way" << std::endl;
    std::cout << "Pos: " << flux.pos.X() << " " << flux.pos.Y() << " " << flux.pos.Z() << std::endl;
    std::cout << "A: " << A.X() << " " << A.Y() << " " << A.Z() << std::endl;
    std::cout << "P: " << flux.mom.Vect().Unit().X() << " " << flux.mom.Vect().Unit().Y() << " " << flux.mom.Vect().Unit().Z() << std::endl;
    return false;
  }

  // set things
  weight = 1. / i;
  flux.mom = higgs_mom;

  if ((flux.pos.Vect() - A).Mag() < (flux.pos.Vect() - B).Mag()) {
    intersection = {A, B}; // A is entry, B is exit
  }
  else {
    intersection = {B, A}; // reversed
  }

  return true;
}

DEFINE_ART_CLASS_TOOL(ReThrowRayTraceBox)

} // namespace ldm
} // namespace evgen
