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

// std includes
#include <string>
#include <iostream>
#include <memory>

//------------------------------------------------------------------------------------------------------------------------------------------
// implementation follows

namespace evgen {
namespace ldm {


// Helper struct retruned by some internal functions
struct RayWeightInfo {
  std::vector<std::vector<TVector3>> allIntersections;
  std::vector<TLorentzVector> allPrtlMom;
  double weight;
  bool pass;
};

/**
 *  @brief  MixedWeightRayTraceBox class definiton
 */
class MixedWeightRayTraceBox : public IRayTrace
{
public:
    /**
     *  @brief  Constructor
     */
    MixedWeightRayTraceBox(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~MixedWeightRayTraceBox();

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
  double fMaxWeightFudge;

  unsigned fNThrow;
  unsigned fNSuccess;
  bool fFixNSuccess;

  double fMaxWeight;

  void CalculateMaxWeight();
  std::pair<double, double> DeltaPhi(TVector3 origin, TRotation &R);
  TLorentzVector ThrowMeVPrtlMomentum(const MeVPrtlFlux &flux, TRotation &RInv, double phi);
  RayWeightInfo ThrowFixedSuccess(const MeVPrtlFlux &flux, TRotation &RInv, double phi);
  RayWeightInfo ThrowFixedThrows(const MeVPrtlFlux &flux, TRotation &RInv, double phi);

};

// Helpers
//
// Minimum Variance Unbiased Estimator for the Negative Binomial Distribution
double QEstimator(unsigned nsuccess, unsigned nfail, unsigned r) {
  if (nsuccess < r) return 0.;

  return ((double)(r - 1) / (r+nfail-1));
}

MixedWeightRayTraceBox::MixedWeightRayTraceBox(fhicl::ParameterSet const &pset):
  IMeVPrtlStage("MixedWeightRayTraceBox") 
{
    this->configure(pset);
}

//------------------------------------------------------------------------------------------------------------------------------------------

MixedWeightRayTraceBox::~MixedWeightRayTraceBox()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
void MixedWeightRayTraceBox::configure(fhicl::ParameterSet const &pset)
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

  fMaxWeightFudge = pset.get<double>("MaxWeightFudge", 1.);

  CalculateMaxWeight();

  fNThrow = pset.get<unsigned>("NThrow");
  fFixNSuccess = pset.get<bool>("FixNSuccess");
  fNSuccess = pset.get<unsigned>("NSuccess", 2);
}
  
void MixedWeightRayTraceBox::CalculateMaxWeight() {
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
  //
  // Cap at 100% intersection probability
  fMaxWeight = std::min(fMaxWeightFudge * weight * fReferenceLabSolidAngle / (4*M_PI), 1.);
}

TLorentzVector MixedWeightRayTraceBox::ThrowMeVPrtlMomentum(const MeVPrtlFlux &flux, TRotation &RInv, double phi) {
  // Pick a direction vector with the parent-direction along the z-axis
  double costh = GetRandom()*2 - 1.;
  double sinth = sqrt(1 - costh*costh);

  TVector3 dir(cos(phi)*sinth, sin(phi)*sinth, costh);

  // Rotate to the lab direction
  dir = RInv * dir;

  // make the mevprtl momentum this
  TLorentzVector mevprtl_mom = flux.mom;
  mevprtl_mom.Boost(-flux.kmom.BoostVector());

  mevprtl_mom = TLorentzVector(mevprtl_mom.P() * dir, mevprtl_mom.E());

  // boost back
  mevprtl_mom.Boost(flux.kmom.BoostVector());
  return mevprtl_mom;
}

std::pair<double, double> MixedWeightRayTraceBox::DeltaPhi(TVector3 origin, TRotation &R) {
  // Iterate over all corners and compute the limits of Phi
  double phihi = -M_PI;
  double philo = M_PI;

  for (unsigned i = 0; i < 8; i++) {
    double X = (i & 1) ? fBox.MaxX() : fBox.MinX();
    double Y = (i & 2) ? fBox.MaxY() : fBox.MinY();
    double Z = (i & 4) ? fBox.MaxZ() : fBox.MinZ();

    TVector3 corner(X, Y, Z);

    double phi = (R * (corner - origin)).Phi();
    if (phi > phihi) phihi = phi;
    if (phi < philo) philo = phi;
  }

  // Now figure out if the bounds should go philo -> phihi or phihi -> philo
  //
  // Whichever has the smaller length is correct
  bool order_is_up = (phihi - philo) < (M_PI - phihi) + (philo + M_PI);

  // If we should go phihi -> philo switch them and re-order them around
  if (!order_is_up) {
    std::cout << "FLIP PHILIM! " << philo << " " << phihi << std::endl;
    double temp = phihi;
    phihi = philo + 2*M_PI;
    philo = temp;
  }

  return {philo, phihi};
}

RayWeightInfo MixedWeightRayTraceBox::ThrowFixedThrows(const MeVPrtlFlux &flux, TRotation &RInv, double phi) {
  RayWeightInfo ret;
  unsigned ithrow = 0;
  unsigned nfail = 0;
  unsigned nsuccess = 0;

  while (ithrow < fNThrow) {
    ithrow ++;

    TLorentzVector mevprtl_mom = ThrowMeVPrtlMomentum(flux, RInv, phi);
    std::vector<TVector3> box_intersections = fBox.GetIntersections(flux.pos.Vect(), mevprtl_mom.Vect().Unit());

    // Does this ray intersect the box?
    if (box_intersections.size() != 2) {
      nfail ++;
      continue;
    }

    // if the ray points the wrong way, it doesn't intersect
    TVector3 A = box_intersections[0];
    if (mevprtl_mom.Vect().Unit().Dot((A - flux.pos.Vect()).Unit()) < 0.) {
      nfail ++;
      continue;
    }

    // if we're here, we have a valid ray
    ret.allIntersections.push_back(box_intersections);
    ret.allPrtlMom.push_back(mevprtl_mom);
    nsuccess ++;
  }

  std::cout << "NTHROW: " << fNThrow << std::endl;
  std::cout << "NSUCCESS: " << nsuccess << std::endl;
  std::cout << "P: " << (((double)nsuccess) / fNThrow) << std::endl;

  ret.weight = ((double)nsuccess) / fNThrow;
  ret.pass = nsuccess > 0;

  return ret;

}

RayWeightInfo MixedWeightRayTraceBox::ThrowFixedSuccess(const MeVPrtlFlux &flux, TRotation &RInv, double phi) {
  RayWeightInfo ret;
  unsigned ithrow = 0;
  unsigned nfail = 0;
  unsigned nsuccess = 0;
  while (nsuccess < fNSuccess && ithrow < fNThrow) {
    ithrow ++;

    TLorentzVector mevprtl_mom = ThrowMeVPrtlMomentum(flux, RInv, phi);
    std::vector<TVector3> box_intersections = fBox.GetIntersections(flux.pos.Vect(), mevprtl_mom.Vect().Unit());

    // Does this ray intersect the box?
    if (box_intersections.size() != 2) {
      nfail ++;
      continue;
    }

    // if the ray points the wrong way, it doesn't intersect
    TVector3 A = box_intersections[0];
    if (mevprtl_mom.Vect().Unit().Dot((A - flux.pos.Vect()).Unit()) < 0.) {
      nfail ++;
      continue;
    }

    // if we're here, we have a valid ray
    ret.allIntersections.push_back(box_intersections);
    ret.allPrtlMom.push_back(mevprtl_mom);
    nsuccess ++;
  }

  std::cout << "NFAIL: " << nfail << std::endl;
  std::cout << "NSUCCESS: " << nsuccess << std::endl;
  std::cout << "Q: " << QEstimator(nsuccess, nfail, fNSuccess) << std::endl;

  if (nsuccess == fNSuccess) {
    ret.weight = QEstimator(nsuccess, nfail, fNSuccess);
    ret.pass = true;
  }
  else {
    ret.weight = 0.;
    ret.pass = false;
  }

  return ret;
}

bool MixedWeightRayTraceBox::IntersectDetector(MeVPrtlFlux &flux, std::array<TVector3, 2> &intersection, double &weight) {
  // We want to pick an orientation in the lab frame (th, phi)
  // that intersects with the detector. 
  //
  // In the lab frame, the angular distribution is not uniform. In the
  // parent-rest-rame (th', phi'), it is. Picking directions in the 
  // parent-rest-frame until a direction which hits the detector is very inefficient.
  //
  // We can draw from the lab frame and weight using d(th, phi) / d(th', phi').
  // (This is what is done for neutrino simulation). However, in the massive
  // case there is a (integrable) singularity which leads to an unbounded max weight.
  // This is undesirable for the Monte-Carlo and makes application of a "accept-reject"
  // algorithm to de-weight impossible. 
  //
  // So here we apply a mixed algorithm. We pick (phi, th'). I.e. we pick phi in the lab
  // frame (forced to hit the detector) and th' in the parent-rest-frame. This reduces 
  // the MC to one dimension instead of two, and so should be much more efficient. And,
  // to weight we just need dphi/dphi' = 1.

  // Setup the axes so that the parent Momentum is along the z axis
  TVector3 parent_dir = flux.kmom.Vect().Unit();
  TVector3 zdir(0, 0, 1);
  TVector3 perp_parent_dir = zdir.Cross(parent_dir);
  double angle = acos(parent_dir.Z());

  TRotation R;

  // Setup the rotation if we need to
  if (perp_parent_dir.Mag() > 1e-4) {
    R.Rotate(-angle, perp_parent_dir);
  }
  TRotation RInv = R.Inverse();

  std::cout << "PARENT DIR: " << parent_dir.X() << " " << parent_dir.Y() << " " << parent_dir.Z() << std::endl;
  TVector3 parent_dir_rotated = R * parent_dir;
  std::cout << "PARENT DIR ROTATED: " << parent_dir_rotated.X() << " " << parent_dir_rotated.Y() << " " << parent_dir_rotated.Z() << std::endl;

  // Compute the Delta Phi of the detector
  std::pair<double, double> philim = DeltaPhi(flux.pos.Vect(), R);
  double philo = philim.first;
  double phihi = philim.second;

  std::cout << "PHILIM: " << philo << " " << phihi << std::endl;
  std::cout << "DELTAPHI: " << (phihi - philo) << std::endl;

  // Pick a phi
  double phi = GetRandom() * (phihi - philo) + philo;

  std::cout << "THISPHI: " << phi << std::endl;

  RayWeightInfo info = (fFixNSuccess) ? ThrowFixedSuccess(flux, RInv, phi) : ThrowFixedThrows(flux, RInv, phi);

  if (!info.pass) return false;

  unsigned ind = CLHEP::RandFlat::shootInt(fEngine, 0, info.allIntersections.size()-1); // inclusive?
  TLorentzVector mevprtl_mom = info.allPrtlMom[ind];
  std::vector<TVector3> box_intersections = info.allIntersections[ind];

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

  // weight from forcing phi
  double phiweight = (phihi - philo) / (2*M_PI);
  // Total weight
  weight = phiweight * info.weight;

  std::cout << "Final WEIGHT: " << weight << std::endl;

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

DEFINE_ART_CLASS_TOOL(MixedWeightRayTraceBox)

} // namespace ldm
} // namespace evgen
