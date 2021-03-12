#ifndef PlaneTransform_HH
#define PlaneTransform_HH

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "canvas/Persistency/Common/FindManyP.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackTrajectory.h"
#include "lardataobj/RecoBase/Trajectory.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larreco/Calorimetry/CalorimetryAlg.h"
#include "sbnobj/Common/Reco/VertexHit.h"
#include "sbnobj/Common/Reco/Stub.h"

#include <memory>

namespace sbn {

// Code for transforming plane coordinates into 3D coordinates
//
// Copied liberally from the Pandora LArRotationalTransformationPlugin
class PlaneTransform {
public:
  PlaneTransform(fhicl::ParameterSet const& p);
  double UVtoW(const double u, const double v) const;
  double VWtoU(const double v, const double w) const;
  double WUtoV(const double w, const double u) const;
  
  double UVtoY(const double u, const double v) const;
  double UVtoZ(const double u, const double v) const;
  double UWtoY(const double u, const double w) const;
  double UWtoZ(const double u, const double w) const;
  double VWtoY(const double v, const double w) const;
  double VWtoZ(const double v, const double w) const;
  
  double YZtoU(const double y, const double z) const;
  double YZtoV(const double y, const double z) const;
  double YZtoW(const double y, const double z) const;

  int ViewtoUVW(geo::View_t v) const;
  double TwoPlaneToY(geo::View_t v1, double w1, geo::View_t v2, double w2) const;
  double TwoPlaneToZ(geo::View_t v1, double w1, geo::View_t v2, double w2) const;
  double YZtoPlane(geo::View_t v, double y, double z) const;
  double WireCoordinate(const geo::GeometryCore *geo, const geo::WireID &w) const;

  void GetMinChiSquaredYZ(const double u, const double v, const double w, const double sigmaU, const double sigmaV, const double sigmaW,
    double &y, double &z, double &chiSquared) const;
  void GetMinChiSquaredYZ(const double u, const double v, const double w, const double sigmaU, const double sigmaV, const double sigmaW,
    const double uFit, const double vFit, const double wFit, const double sigmaFit, double &y, double &z, double &chiSquared) const;
  

private:

  double    m_thetaU;                 ///< inclination of U wires (radians)
  double    m_thetaV;                 ///< inclination of V wires (radians)
  double    m_thetaW;                 ///< inclination of W wires (radians)
  
  double    m_sinU;                   ///< sin(thetaU)
  double    m_sinV;                   ///< sin(thetaV)
  double    m_sinW;                   ///< sin(thetaW)
  double    m_cosU;                   ///< cos(thetaU)
  double    m_cosV;                   ///< cos(thetaV)
  double    m_cosW;                   ///< cos(thetaW)
  double    m_sinVminusU;             ///< sin(thetaV - thetaU)
  double    m_sinWminusV;             ///< sin(thetaW - thetaV)
  double    m_sinUminusW;             ///< sin(thetaU - thetaW)
  
  double    m_maxAngularDiscrepancyU; ///< Maximum allowed difference between u wire angles between LArTPCs
  double    m_maxAngularDiscrepancyV; ///< Maximum allowed difference between v wire angles between LArTPCs
  double    m_maxAngularDiscrepancyW; ///< Maximum allowed difference between w wire angles between LArTPCs
  double    m_maxSigmaDiscrepancy;    ///< Maximum allowed difference between like wire sigma values between LArTPCs

  std::vector<geo::View_t> m_vieworder;

};

} // end namespace sbn
#endif
