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
#include "../Products/MeVPrtlFlux.h"
#include "sbncode/EventGenerator/MeVPrtl/Tools/Constants.h"

// LArSoft includes
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/CoreUtils/ServiceUtil.h"

// for testing
#include "dk2nu/tree/calcLocationWeights.h"
#include "dk2nu/tree/dkmeta.h"
#include "dk2nu/tree/dk2nu.h"

// std includes
#include <string>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace evgen {
namespace ldm {

/**
 *  @brief  WeightedRayTraceBox class definiton
 */
class WeightedRayTraceBox : public IRayTrace
{
public:
    /**
     *  @brief  Constructor
     */
    WeightedRayTraceBox(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~WeightedRayTraceBox();

    void configure(const fhicl::ParameterSet&) override;

    bool IntersectDetector(MeVPrtlFlux &flux, std::array<TVector3, 2> &intersection, double &weight) override;

    // always thrown at least once
    float MaxWeight() override { 
      return fMaxWeight;
    }

private:
  geo::BoxBoundedGeo fBox;
  double fReferenceLabSolidAngle;
  double fReferencePrtlMass;
  int fReferenceScndPDG;
  double fReferenceKaonEnergy;

  void CalculateMaxWeight();

  double fMaxWeight;

  double SolidAngle(TVector3 loc);
  TVector3 RandomIntersectionPoint(TVector3 loc);
};

WeightedRayTraceBox::WeightedRayTraceBox(fhicl::ParameterSet const &pset):
  IMeVPrtlStage("WeightedRayTraceBox") 
{
    this->configure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

WeightedRayTraceBox::~WeightedRayTraceBox()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
void WeightedRayTraceBox::configure(fhicl::ParameterSet const &pset)
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

  std::cout << "Detector Box." << std::endl;
  std::cout << "X " << fBox.MinX() << " " << fBox.MaxX() << std::endl;
  std::cout << "Y " << fBox.MinY() << " " << fBox.MaxY() << std::endl;
  std::cout << "Z " << fBox.MinZ() << " " << fBox.MaxZ() << std::endl;

  fReferenceLabSolidAngle = pset.get<double>("ReferenceLabSolidAngle");
  fReferencePrtlMass = pset.get<double>("ReferencePrtlMass");
  fReferenceScndPDG = pset.get<int>("ReferenceScndPDG");
  fReferenceKaonEnergy = pset.get<double>("ReferenceKaonEnergy");

  CalculateMaxWeight();

  // TEST: compare this weight calc to calcENuWeight
  // uncomment to run
  /*
  double p = twobody_momentum(kaonp_mass, pionp_mass, 0.);

  {
    double beta = sqrt(fReferenceKaonEnergy*fReferenceKaonEnergy - kaonp_mass*kaonp_mass) / fReferenceKaonEnergy;
    // Doen't affect weight
    double rand = 1.;

    // For the maxweight, make the kaon and daughter direction aligned
    TVector3 dir(0., 0, 1.);
    TVector3 boost(0., 0., beta);
    double lab_frame_p, costh_rest, weight;
    int err = calcPrtlRayWgt(p, 0., dir, boost, rand, lab_frame_p, costh_rest, weight);
    (void) err;

    std::cout << "REST FRAME P: " << p << std::endl;
    std::cout << "THIS WEIGHT: " << weight << std::endl;
    std::cout << "THIS P: " << lab_frame_p << std::endl;
  }

  {
    TVector3 start(0, 0, 0);
    TVector3 end(0, 0, 100e2);

    bsim::Decay decay;
    decay.ptype = 321;
    decay.pdpx = 0.;
    decay.pdpy = 0.;
    decay.pdpz = sqrt(fReferenceKaonEnergy*fReferenceKaonEnergy - kaonp_mass*kaonp_mass);
    decay.necm = p;
    decay.vx = start.X();
    decay.vy = start.Y();
    decay.vz = start.Z();
    decay.ndecay = 0;

    double weight, enu;
    int err = bsim::calcEnuWgt(decay, end, enu, weight);
    (void) err;

    std::cout << "BSIM WEIGHT: " << weight << std::endl;
    std::cout << "BSIM P: " << enu << std::endl;
  }
  */

}
  
void WeightedRayTraceBox::CalculateMaxWeight() {
  double secondary_mass = secPDG2Mass(fReferenceScndPDG);

  double p = twobody_momentum(kaonp_mass, secondary_mass, fReferencePrtlMass);
  double beta = sqrt(fReferenceKaonEnergy*fReferenceKaonEnergy - kaonp_mass*kaonp_mass) / fReferenceKaonEnergy;

  // Doen't affect weight
  double rand = 1.;

  // For the maxweight, make the kaon and daughter direction aligned
  TVector3 dir(1, 0, 0);
  TVector3 boost(beta, 0., 0.);

  double lab_frame_p, costh_rest, weight;
  int err = calcPrtlRayWgt(p, fReferencePrtlMass, dir, boost, rand, lab_frame_p, costh_rest, weight);

  if (err) { // Shouldn't happen
    throw cet::exception("Weighted Ray Trace: Bad max weight calculation.");
  }

  // Normalize the weight
  fMaxWeight = weight * fReferenceLabSolidAngle / (4*M_PI);
}

TVector3 WeightedRayTraceBox::RandomIntersectionPoint(TVector3 loc) {
  // In each dimension, figure out which of the face the location would intesect
  bool intersect_MinX = loc.X() < fBox.MinX();
  bool intersect_MinY = loc.Y() < fBox.MinY();
  bool intersect_MinZ = loc.Z() < fBox.MinZ();

  bool intersect_MaxX = loc.X() > fBox.MaxX();
  bool intersect_MaxY = loc.Y() > fBox.MaxY();
  bool intersect_MaxZ = loc.Z() > fBox.MaxZ();
  
  // Solid angle of each surface
  double solid_angleX = 0.;
  double solid_angleY = 0.;
  double solid_angleZ = 0.;

  // Compute the contribution to the solid angle of each face
  if (intersect_MinX) {
    double area = (fBox.MaxY() - fBox.MinY()) * (fBox.MaxZ() - fBox.MinZ());
    TVector3 center(fBox.MinX(), fBox.CenterY(), fBox.CenterZ());

    solid_angleX = area * abs((center-loc).Unit().X()) / (center-loc).Mag2();
  }
  else if (intersect_MaxX) {
    double area = (fBox.MaxY() - fBox.MinY()) * (fBox.MaxZ() - fBox.MinZ());
    TVector3 center(fBox.MaxX(), fBox.CenterY(), fBox.CenterZ());

    solid_angleX = area * abs((center-loc).Unit().X()) / (center-loc).Mag2();
  }

  if (intersect_MinY) {
    double area = (fBox.MaxX() - fBox.MinX()) * (fBox.MaxZ() - fBox.MinZ());
    TVector3 center(fBox.CenterX(), fBox.MinY(), fBox.CenterZ());

    solid_angleY = area * abs((center-loc).Unit().Y()) / (center-loc).Mag2();
  }
  else if (intersect_MaxY) {
    double area = (fBox.MaxX() - fBox.MinX()) * (fBox.MaxZ() - fBox.MinZ());
    TVector3 center(fBox.CenterX(), fBox.MaxY(), fBox.CenterZ());

    solid_angleY = area * abs((center-loc).Unit().Y()) / (center-loc).Mag2();
  }

  if (intersect_MinZ) {
    double area = (fBox.MaxY() - fBox.MinY()) * (fBox.MaxX() - fBox.MinX());
    TVector3 center(fBox.CenterX(), fBox.CenterY(), fBox.MinZ());

    solid_angleZ = area * abs((center-loc).Unit().Z()) / (center-loc).Mag2();
  }
  else if (intersect_MaxZ) {
    double area = (fBox.MaxY() - fBox.MinY()) * (fBox.MaxX() - fBox.MinX());
    TVector3 center(fBox.CenterX(), fBox.CenterY(), fBox.MaxZ());

    solid_angleZ = area * abs((center-loc).Unit().Z()) / (center-loc).Mag2();
  }

  // Pick which face to intersect
  double rand = GetRandom();

  bool select_X = rand < solid_angleX / (solid_angleX + solid_angleY + solid_angleZ);
  bool select_Y = !select_X && (rand < (solid_angleX + solid_angleY) / (solid_angleX + solid_angleY + solid_angleZ));
  bool select_Z = !select_X && !select_Y;

  double X, Y, Z;
  if (select_X) {
    X = (intersect_MinX) ? fBox.MinX() : fBox.MaxX(); 
    Y = (fBox.MaxY() - fBox.MinY()) * GetRandom() + fBox.MinY();
    Z = (fBox.MaxZ() - fBox.MinZ()) * GetRandom() + fBox.MinZ();
  }
  if (select_Y) {
    X = (fBox.MaxX() - fBox.MinX()) * GetRandom() + fBox.MinX();
    Y = (intersect_MinY) ? fBox.MinY() : fBox.MaxY(); 
    Z = (fBox.MaxZ() - fBox.MinZ()) * GetRandom() + fBox.MinZ();
  }
  if (select_Z) {
    X = (fBox.MaxX() - fBox.MinX()) * GetRandom() + fBox.MinX();
    Y = (fBox.MaxY() - fBox.MinY()) * GetRandom() + fBox.MinY();
    Z = (intersect_MinZ) ? fBox.MinZ() : fBox.MaxZ(); 
  }

  TVector3 ret(X, Y, Z);
  return ret;
}

double WeightedRayTraceBox::SolidAngle(TVector3 loc) {
  // In each dimension, figure out which of the face the location would intesect
  bool intersect_MinX = loc.X() < fBox.MinX();
  bool intersect_MinY = loc.Y() < fBox.MinY();
  bool intersect_MinZ = loc.Z() < fBox.MinZ();

  bool intersect_MaxX = loc.X() > fBox.MaxX();
  bool intersect_MaxY = loc.Y() > fBox.MaxY();
  bool intersect_MaxZ = loc.Z() > fBox.MaxZ();
  
  double solid_angle = 0.;

  // Compute the contribution to the solid angle of each face
  if (intersect_MinX) {
    double area = (fBox.MaxY() - fBox.MinY()) * (fBox.MaxZ() - fBox.MinZ());
    TVector3 center(fBox.MinX(), fBox.CenterY(), fBox.CenterZ());

    solid_angle += area * abs((center-loc).Unit().X()) / (center-loc).Mag2();
  }
  if (intersect_MinY) {
    double area = (fBox.MaxX() - fBox.MinX()) * (fBox.MaxZ() - fBox.MinZ());
    TVector3 center(fBox.CenterX(), fBox.MinY(), fBox.CenterZ());

    solid_angle += area * abs((center-loc).Unit().Y()) / (center-loc).Mag2();
  }
  if (intersect_MinZ) {
    double area = (fBox.MaxY() - fBox.MinY()) * (fBox.MaxX() - fBox.MinX());
    TVector3 center(fBox.CenterX(), fBox.CenterY(), fBox.MinZ());

    solid_angle += area * abs((center-loc).Unit().Z()) / (center-loc).Mag2();
  }
  
  if (intersect_MaxX) {
    double area = (fBox.MaxY() - fBox.MinY()) * (fBox.MaxZ() - fBox.MinZ());
    TVector3 center(fBox.MaxX(), fBox.CenterY(), fBox.CenterZ());

    solid_angle += area * abs((center-loc).Unit().X()) / (center-loc).Mag2();
  }
  if (intersect_MaxY) {
    double area = (fBox.MaxX() - fBox.MinX()) * (fBox.MaxZ() - fBox.MinZ());
    TVector3 center(fBox.CenterX(), fBox.MaxY(), fBox.CenterZ());

    solid_angle += area * abs((center-loc).Unit().Y()) / (center-loc).Mag2();
  }
  if (intersect_MaxZ) {
    double area = (fBox.MaxY() - fBox.MinY()) * (fBox.MaxX() - fBox.MinX());
    TVector3 center(fBox.CenterX(), fBox.CenterY(), fBox.MaxZ());

    solid_angle += area * abs((center-loc).Unit().Z()) / (center-loc).Mag2();
  }

  return solid_angle;

}
    
bool WeightedRayTraceBox::IntersectDetector(MeVPrtlFlux &flux, std::array<TVector3, 2> &intersection, double &weight) {
  // Randomly pick a location in the detector to send this particle to
  // TVector3 detloc = RandomIntersectionPoint(flux.pos.Vect());

  // For a box far away from the origin, picking a point uniformly in the volume of the box
  // is the same as picking a point uniformly on the solid angle about that origin
  double detX = (fBox.MaxX() - fBox.MinX()) * GetRandom() + fBox.MinX();
  double detY = (fBox.MaxY() - fBox.MinY()) * GetRandom() + fBox.MinY();
  double detZ = (fBox.MaxZ() - fBox.MinZ()) * GetRandom() + fBox.MinZ();
  TVector3 detloc(detX, detY, detZ);

  // calculate the weight and kinematics of sending the flux here
  TLorentzVector flux_mom_rest = flux.mom;
  flux_mom_rest.Boost(-flux.kmom.BoostVector());
  double rest_frame_p = flux_mom_rest.Vect().Mag();
  double p_lab, costh_rest;
  int err = calcPrtlRayWgt(rest_frame_p, flux.mom.M(), detloc - flux.pos.Vect(), flux.kmom.BoostVector(), GetRandom(), 
      p_lab, costh_rest, weight);

  // Kinematics don't work
  if (err) return false;

  // Adjust the input flux to the computed kinematics
  flux.mom.SetVectM(p_lab * (detloc - flux.pos.Vect()).Unit(), flux.mom.M());

  // And correct all the other variables
  flux.mom_beamcoord = flux.mom;
  flux.mom_beamcoord.Boost(-flux.kmom.BoostVector());
  flux.mom_beamcoord.Boost(flux.kmom_beamcoord.BoostVector());

  flux.sec = flux.kmom - flux.mom;
  flux.sec_beamcoord = flux.kmom_beamcoord - flux.mom_beamcoord;

  // Compute the intersections for the selected point
  std::vector<TVector3> box_intersections = fBox.GetIntersections(flux.pos.Vect(), flux.mom.Vect().Unit());
  TVector3 A = box_intersections.at(0);
  TVector3 B = box_intersections.at(1);
  if ((flux.pos.Vect() - A).Mag() < (flux.pos.Vect() - B).Mag()) {
    intersection = {A, B}; // A is entry, B is exit
  }
  else {
    intersection = {B, A}; // reversed
  }

  // Turn the weight into an event weight
  weight *= SolidAngle(flux.pos.Vect()) / (4*M_PI);

  std::cout << "From: " << flux.pos.X() << " " << flux.pos.Y() << " " << flux.pos.Z() << std::endl;
  std::cout << "Solid Angle ratio is: " << (SolidAngle(flux.pos.Vect()) / (4*M_PI)) << std::endl;

  std::cout << "Kaon 4P: " << flux.kmom.E() << " " << flux.kmom.Px() << " " << flux.kmom.Py() << " " << flux.kmom.Pz() << std::endl;
  std::cout << "Selected Prtl 4P: " << flux.mom.E() << " " << flux.mom.Px() << " " << flux.mom.Py() << " " << flux.mom.Pz() << std::endl;
  std::cout << "Selected Scdy 4P: " << flux.sec.E() << " " << flux.sec.Px() << " " << flux.sec.Py() << " " << flux.sec.Pz() << std::endl;

  // check the costh calc
  TLorentzVector flux_rest = flux.mom;
  flux_rest.Boost(-flux.kmom.BoostVector());
  std::cout << "Prtl lab costh: " << costh_rest << " " << flux.kmom.Vect().Unit().Dot(flux_rest.Vect().Unit()) << std::endl;

  return true;
}

DEFINE_ART_CLASS_TOOL(WeightedRayTraceBox)

} // namespace ldm
} // namespace evgen
