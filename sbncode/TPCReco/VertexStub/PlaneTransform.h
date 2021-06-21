#ifndef PlaneTransform_HH
#define PlaneTransform_HH

#include "larcorealg/Geometry/GeometryCore.h"

#include <vector>

namespace sbn {

// Code for transforming plane coordinates into 3D coordinates
//
// Copied liberally from the Pandora LArRotationalTransformationPlugin
class PlaneTransform {
public:
  PlaneTransform(fhicl::ParameterSet const& p, const geo::GeometryCore *geo);
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

  int PlanetoUVW(const geo::PlaneID &p) const;
  double TwoPlaneToY(const geo::PlaneID &p1, double w1, const geo::PlaneID &p2, double w2) const;
  double TwoPlaneToZ(const geo::PlaneID &p1, double w1, const geo::PlaneID &p2, double w2) const;
  double YZtoPlane(const geo::PlaneID &p, double y, double z) const;
  double WireCoordinate(const geo::WireID &w) const;

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

  const     geo::GeometryCore *fGeo;   ///< Handle to the geometry
};

} // end namespace sbn
#endif
