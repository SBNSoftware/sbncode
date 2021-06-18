#include "PlaneTransform.h"

sbn::PlaneTransform::PlaneTransform(fhicl::ParameterSet const& p):
  m_thetaU(p.get<double>("WireAngleU")),
  m_thetaV(p.get<double>("WireAngleV")),
  m_thetaW(p.get<double>("WireAngleW")),
  m_sinU(sin(m_thetaU)),
  m_sinV(sin(m_thetaV)),
  m_sinW(sin(m_thetaW)),
  m_cosU(cos(m_thetaU)),
  m_cosV(cos(m_thetaV)),
  m_cosW(cos(m_thetaW)),
  m_sinVminusU(sin(m_thetaV - m_thetaU)),
  m_sinWminusV(sin(m_thetaW - m_thetaV)),
  m_sinUminusW(sin(m_thetaU - m_thetaW))
{
  if (abs(m_sinVminusU) < std::numeric_limits<double>::epsilon() ||
      abs(m_sinWminusV) < std::numeric_limits<double>::epsilon() ||
      abs(m_sinUminusW) < std::numeric_limits<double>::epsilon()) {
    throw cet::exception("PlaneTransform::PlaneTransform: Does not support provided TPC configuration.");
  }

  std::vector<unsigned> ViewOrderConfig = p.get<std::vector<unsigned>>("ViewOrder");
  for (unsigned v: ViewOrderConfig) m_vieworder.push_back((geo::View_t)v);

}

double sbn::PlaneTransform::UVtoW(const double u, const double v) const
{
  return (-1. * (u * m_sinWminusV + v * m_sinUminusW) / m_sinVminusU);
}


double sbn::PlaneTransform::VWtoU(const double v, const double w) const
{
  return (-1. * (v * m_sinUminusW + w * m_sinVminusU) / m_sinWminusV);
}

double sbn::PlaneTransform::WUtoV(const double w, const double u) const
{
  return (-1. * (u * m_sinWminusV + w * m_sinVminusU) / m_sinUminusW);
}


double sbn::PlaneTransform::UVtoY(const double u, const double v) const
{
  return ((u * m_cosV - v * m_cosU) / m_sinVminusU);
}

double sbn::PlaneTransform::UVtoZ(const double u, const double v) const
{
  return ((u * m_sinV - v * m_sinU) / m_sinVminusU);
}

double sbn::PlaneTransform::UWtoY(const double u, const double w) const
{
  return ((w * m_cosU - u * m_cosW) / m_sinUminusW);
}

double sbn::PlaneTransform::UWtoZ(const double u, const double w) const
{
  return ((w * m_sinU - u * m_sinW) / m_sinUminusW);
}

double sbn::PlaneTransform::VWtoY(const double v, const double w) const
{
  return ((v * m_cosW - w * m_cosV) / m_sinWminusV);
}

double sbn::PlaneTransform::VWtoZ(const double v, const double w) const
{
  return ((v * m_sinW - w * m_sinV) / m_sinWminusV);
}

double sbn::PlaneTransform::YZtoU(const double y, const double z) const
{
  return (z * m_cosU - y * m_sinU);
}

double sbn::PlaneTransform::YZtoV(const double y, const double z) const
{
  return (z * m_cosV - y * m_sinV);
}

double sbn::PlaneTransform::YZtoW(const double y, const double z) const
{
  return (z * m_cosW - y * m_sinW);
}

int sbn::PlaneTransform::ViewtoUVW(geo::View_t v) const {
  for (unsigned i = 0; i < m_vieworder.size(); i++) {
    if (m_vieworder[i] == v) return i;
  }
  return -1;
}

double sbn::PlaneTransform::TwoPlaneToY(geo::View_t v1, double w1, geo::View_t v2, double w2) const {
  unsigned uvw1 = ViewtoUVW(v1);
  unsigned uvw2 = ViewtoUVW(v2);

  if (uvw1 == 0 && uvw2 == 1) {
    return UVtoY(w1, w2);
  }
  else if (uvw1 == 1 && uvw2 == 0) {
    return UVtoY(w2, w1);
  }
  else if (uvw1 == 1 && uvw2 == 2) {
    return VWtoY(w1, w2);
  }
  else if (uvw1 == 2 && uvw2 == 1) {
    return VWtoY(w2, w1);
  }
  else if (uvw1 == 0 && uvw2 == 2) {
    return UWtoY(w1, w2);
  }
  else if (uvw1 == 2 && uvw2 == 0) {
    return UWtoY(w2, w1);
  }
  else {} // BAD!

  return -100000.;
}

double sbn::PlaneTransform::TwoPlaneToZ(geo::View_t v1, double w1, geo::View_t v2, double w2) const {
  unsigned uvw1 = ViewtoUVW(v1);
  unsigned uvw2 = ViewtoUVW(v2);

  if (uvw1 == 0 && uvw2 == 1) {
    return UVtoZ(w1, w2);
  }
  else if (uvw1 == 1 && uvw2 == 0) {
    return UVtoZ(w2, w1);
  }
  else if (uvw1 == 1 && uvw2 == 2) {
    return VWtoZ(w1, w2);
  }
  else if (uvw1 == 2 && uvw2 == 1) {
    return VWtoZ(w2, w1);
  }
  else if (uvw1 == 0 && uvw2 == 2) {
    return UWtoZ(w1, w2);
  }
  else if (uvw1 == 2 && uvw2 == 0) {
    return UWtoZ(w2, w1);
  }
  else {} // BAD!

  return -100000.;
}

double sbn::PlaneTransform::YZtoPlane(geo::View_t v, double y, double z) const {
  unsigned uvw = ViewtoUVW(v);

  if (uvw == 0) return YZtoU(y, z);
  else if (uvw == 1) return YZtoV(y, z);
  else if (uvw == 2) return YZtoW(y, z);
  else {} // BAD!

  return -100000.;
}

double sbn::PlaneTransform::WireCoordinate(const geo::GeometryCore *geo, const geo::WireID &w) const {
  geo::View_t v = geo->View(w);
  TVector3 coord = geo->Wire(w).GetCenter();

  return YZtoPlane(v, coord.Y(), coord.Z());
}

