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
#include "sbncode/EventGenerator/MeVPrtl/Tools/Constants.h"

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

    bool IntersectDetector(MeVPrtlFlux &flux, std::array<TVector3, 2> &intersection, double &weight) override;

    TLorentzVector ThrowMeVPrtlMomentum(const MeVPrtlFlux &flux);

    // always thrown at least once
    float MaxWeight() override { 
      return fMaxWeight;
    }

private:
  geo::BoxBoundedGeo fBox;
  unsigned fNThrows;

  double fReferenceLabSolidAngle;
  double fReferencePrtlMass;
  int fReferenceScndPDG;
  double fReferenceKaonEnergy;

  void CalculateMaxWeight();

  double fMaxWeight;
};

ReThrowRayTraceBox::ReThrowRayTraceBox(fhicl::ParameterSet const &pset):
  IMeVPrtlStage("ReThrowRayTraceBox") 
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
  if (pset.has_key("Box")) {
    std::array<double, 6> box_config = pset.get<std::array<double, 6>>("Box");
    // xmin, xmax, ymin, ymax, zmin, zmax
    fBox = geo::BoxBoundedGeo(box_config[0], box_config[1], box_config[2], box_config[3], box_config[4], box_config[5]);
  }
  else {
    const geo::GeometryCore *geometry = lar::providerFrom<geo::Geometry>();
    fBox = geometry->DetectorEnclosureBox(pset.get<std::string>("Volume"));
  }

  fNThrows = pset.get<unsigned>("NThrows", 10000);

  std::cout << "Detector Box." << std::endl;
  std::cout << "X " << fBox.MinX() << " " << fBox.MaxX() << std::endl;
  std::cout << "Y " << fBox.MinY() << " " << fBox.MaxY() << std::endl;
  std::cout << "Z " << fBox.MinZ() << " " << fBox.MaxZ() << std::endl;

  fReferenceLabSolidAngle = pset.get<double>("ReferenceLabSolidAngle");
  fReferencePrtlMass = pset.get<double>("ReferencePrtlMass");
  fReferenceScndPDG = pset.get<int>("ReferenceScndPDG");
  fReferenceKaonEnergy = pset.get<double>("ReferenceKaonEnergy");

  CalculateMaxWeight();

  // TEST: compare results of this weight to calcENuWeight

}
  
void ReThrowRayTraceBox::CalculateMaxWeight() {
  double secondary_mass = 0.;
  switch (fReferenceScndPDG) {
    case 11:
    case -11:
      secondary_mass = Constants::Instance().elec_mass;
      break;
    case 13:
    case -13:
      secondary_mass = Constants::Instance().muon_mass;
      break;
    case 211:
    case -211:
      secondary_mass = Constants::Instance().piplus_mass;
      break;
    default:
      std::cerr << "RETHROWRAYTACE -- bad secondary pdg: " << fReferenceScndPDG << std::endl;
      std::cout << "RETHROWRAYTACE -- bad secondary pdg: " << fReferenceScndPDG << std::endl;
      break;
  }

  double p = twobody_momentum(Constants::Instance().kplus_mass, secondary_mass, fReferencePrtlMass);
  double E = sqrt(p*p + fReferencePrtlMass*fReferencePrtlMass);

  double kgamma = fReferenceKaonEnergy / Constants::Instance().kplus_mass;
  double kbeta = sqrt(1 - 1. / (kgamma * kgamma));

  double scale = kgamma * kgamma * (p + kbeta*E) * (p + kbeta*E) / (p*p);

  fMaxWeight = scale * fReferenceLabSolidAngle;
}

TLorentzVector ReThrowRayTraceBox::ThrowMeVPrtlMomentum(const MeVPrtlFlux &flux) {
  // pick a direction for the rest-frame
  TVector3 dir = RandomUnitVector();

  // make the mevprtl momentum this
  // 
  // Move the mevprtl momentum back to the kaon rest frame
  TLorentzVector mevprtl_mom = flux.mom;
  mevprtl_mom.Boost(-flux.kmom.BoostVector());

  mevprtl_mom = TLorentzVector(mevprtl_mom.P() * dir, mevprtl_mom.E());

  // boost back
  mevprtl_mom.Boost(flux.kmom.BoostVector());

  return mevprtl_mom;
  
}
    
bool ReThrowRayTraceBox::IntersectDetector(MeVPrtlFlux &flux, std::array<TVector3, 2> &intersection, double &weight) {
  // try out the mevprtl direction a bunch of times
  std::vector<std::vector<TVector3>> allIntersections;
  std::vector<TLorentzVector> allHMom;

  for (unsigned i = 0; i <= fNThrows; i++) {
    TLorentzVector mevprtl_mom = ThrowMeVPrtlMomentum(flux);
    std::vector<TVector3> box_intersections = fBox.GetIntersections(flux.pos.Vect(), mevprtl_mom.Vect().Unit());

    // Does this ray intersect the box?
    if (box_intersections.size() != 2) {
      continue;
    }

    TVector3 A = box_intersections[0];

    // if the ray points the wrong way, it doesn't intersect
    if (mevprtl_mom.Vect().Unit().Dot((A - flux.pos.Vect()).Unit()) < 0.) {
      std::cout << "RAYTRACE: MeVPrtl points wrong way" << std::endl;
      std::cout << "Pos: " << flux.pos.X() << " " << flux.pos.Y() << " " << flux.pos.Z() << std::endl;
      std::cout << "A: " << A.X() << " " << A.Y() << " " << A.Z() << std::endl;
      std::cout << "P: " << mevprtl_mom.Vect().Unit().X() << " " << mevprtl_mom.Vect().Unit().Y() << " " << mevprtl_mom.Vect().Unit().Z() << std::endl;
      continue;
    }

    // if we're here, we have a valid ray
    allIntersections.push_back(box_intersections);
    allHMom.push_back(mevprtl_mom);
  }

  std::cout << "Prtl intersected (" << allIntersections.size() << " / " << fNThrows << ") times.\n";

  // did we get a hit?
  if (allIntersections.size() == 0) {
    return false;
  }

  unsigned ind = CLHEP::RandFlat::shootInt(fEngine, 0, allIntersections.size()-1); // inclusive?

  TLorentzVector mevprtl_mom = allHMom[ind];
  std::vector<TVector3> box_intersections = allIntersections[ind];

  TVector3 A = box_intersections[0];
  TVector3 B = box_intersections[1];

  // make sure that the flux start lies outside the detector
  if ((flux.pos.Vect() - A).Mag() < (A-B).Mag() && (flux.pos.Vect() - B).Mag() < (A-B).Mag()) {
    throw cet::exception("ReThrowRayTraceBox Exception", "Input mevprtl flux starts inside detector volume: "
        "MeVPrtl start At: (" + std::to_string(flux.pos.X()) + ", " + std::to_string(flux.pos.Y()) + ", " + std::to_string(flux.pos.Z()) + "). "
        "Intersection A At: " + std::to_string(A.X()) + ", " + std::to_string(A.Y()) + ", " + std::to_string(A.Z()) + "). "
        "Intersection B At: " + std::to_string(B.X()) + ", " + std::to_string(B.Y()) + ", " + std::to_string(B.Z()) + ").\n"
    );
  } 

  // set things
  weight = (double)allIntersections.size() / fNThrows;
  flux.mom = mevprtl_mom;
  // transform to beam-coord frame
  flux.mom_beamcoord = mevprtl_mom;
  flux.mom_beamcoord.Boost(-flux.kmom.BoostVector()); // Boost to kaon rest frame
  flux.mom_beamcoord.Boost(flux.kmom_beamcoord.BoostVector()); // And to beam coordinate frame
  // Two boosts make a rotation!

  // Set the secondary 4-momentum
  flux.sec = flux.kmom - flux.mom;
  flux.sec_beamcoord = flux.kmom_beamcoord - flux.mom_beamcoord;

  if ((flux.pos.Vect() - A).Mag() < (flux.pos.Vect() - B).Mag()) {
    intersection = {A, B}; // A is entry, B is exit
  }
  else {
    intersection = {B, A}; // reversed
  }

  std::cout << "Kaon 4P: " << flux.kmom.E() << " " << flux.kmom.Px() << " " << flux.kmom.Py() << " " << flux.kmom.Pz() << std::endl;
  std::cout << "Selected Prtl 4P: " << flux.mom.E() << " " << flux.mom.Px() << " " << flux.mom.Py() << " " << flux.mom.Pz() << std::endl;
  std::cout << "Selected Scdy 4P: " << flux.sec.E() << " " << flux.sec.Px() << " " << flux.sec.Py() << " " << flux.sec.Pz() << std::endl;

  return true;
}

DEFINE_ART_CLASS_TOOL(ReThrowRayTraceBox)

} // namespace ldm
} // namespace evgen
