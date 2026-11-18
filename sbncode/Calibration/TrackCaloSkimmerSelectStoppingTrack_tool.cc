// Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcorealg/Geometry/GeometryCore.h"

#include "ITCSSelectionTool.h"

namespace sbn {

class TrackCaloSkimmerSelectStoppingTrack: public ITCSSelectionTool {
public:

  TrackCaloSkimmerSelectStoppingTrack(const fhicl::ParameterSet &p);
  ~TrackCaloSkimmerSelectStoppingTrack() {}

  bool Select(const TrackInfo &t) override;

private:
  // config
  double fFVInsetMinX;
  double fFVInsetMaxX;
  double fFVInsetMinY;
  double fFVInsetMaxY;
  double fFVInsetMinZ;
  double fFVInsetMaxZ;

  double fMinTimeTickInset;
  double fMaxTimeTickInset;

  double fEndMediandQdxCut;
  unsigned fNumberTimeSamples;
  bool fRequireDownwards;

  // Time info grabbed from Detector Properties
  int fTickMin;
  int fTickMax;

  // Geometry info grabbed from Detector Geometry
  std::vector<std::vector<geo::BoxBoundedGeo>> fTPCVolumes;
  std::vector<geo::BoxBoundedGeo> fActiveVolumes;

  // Fiducial space
  std::vector<geo::BoxBoundedGeo> fFiducialVolumes;
  // Fiducial time
  double fFidTickMin;
  double fFidTickMax;
  double fMediandQdxRRMax;

  bool fCheckFiducialX;
};

TrackCaloSkimmerSelectStoppingTrack::TrackCaloSkimmerSelectStoppingTrack(const fhicl::ParameterSet &p):
  ITCSSelectionTool(p),
  fFVInsetMinX(p.get<double>("FVInsetMinX")),
  fFVInsetMaxX(p.get<double>("FVInsetMaxX")),
  fFVInsetMinY(p.get<double>("FVInsetMinY")),
  fFVInsetMaxY(p.get<double>("FVInsetMaxY")),
  fFVInsetMinZ(p.get<double>("FVInsetMinZ")),
  fFVInsetMaxZ(p.get<double>("FVInsetMaxZ")),
  fMinTimeTickInset(p.get<double>("MinTimeTickInset")),
  fMaxTimeTickInset(p.get<double>("MaxTimeTickInset")),
  fEndMediandQdxCut(p.get<double>("EndMediandQdxCut")),
  fNumberTimeSamples(p.get<unsigned>("NumberTimeSamples")),
  fRequireDownwards(p.get<bool>("RequireDownwards", true)),
  fMediandQdxRRMax(p.get<double>("MediandQdxRRMax", 5.)),
  fCheckFiducialX(p.get<bool>("CheckFiducialX"))
{
  // Get the fiducial volume info
  const geo::GeometryCore *geometry = lar::providerFrom<geo::Geometry>();

  // first the TPC volumes 
  for (auto const &cryoid: geometry->Iterate<geo::CryostatID>()) {
    std::vector<geo::BoxBoundedGeo> this_tpc_volumes;
    for (auto const& TPC : geometry->Iterate<geo::TPCGeo>(cryoid)) {
      this_tpc_volumes.push_back(TPC.ActiveBoundingBox());
    }
     fTPCVolumes.push_back(std::move(this_tpc_volumes));
  }

  // TODO: make configurable? Is this every not 0?
  fTickMin = 0;
  fTickMax = fTickMin + fNumberTimeSamples;

  fFidTickMin = fTickMin + fMinTimeTickInset;
  fFidTickMax = fTickMax - fMaxTimeTickInset;

  // then combine them into active volumes
  for (const std::vector<geo::BoxBoundedGeo> &tpcs: fTPCVolumes) {
    double XMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinX() < rhs.MinX(); })->MinX();
    double YMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinY() < rhs.MinY(); })->MinY();
    double ZMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinZ() < rhs.MinZ(); })->MinZ();

    double XMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxX() < rhs.MaxX(); })->MaxX();
    double YMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxY() < rhs.MaxY(); })->MaxY();
    double ZMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxZ() < rhs.MaxZ(); })->MaxZ();

    fActiveVolumes.emplace_back(XMin, XMax, YMin, YMax, ZMin, ZMax);
  }

  // And take the inset for the fiducial volumes
  for (const geo::BoxBoundedGeo &AV: fActiveVolumes) {
    fFiducialVolumes.emplace_back(AV.MinX() + fFVInsetMinX, AV.MaxX() - fFVInsetMaxX, 
                                  AV.MinY() + fFVInsetMinY, AV.MaxY() - fFVInsetMaxY, 
                                  AV.MinZ() + fFVInsetMinZ, AV.MaxZ() - fFVInsetMaxZ);
  }

}

bool TrackCaloSkimmerSelectStoppingTrack::Select(const TrackInfo &t) {
  bool downwards = (t.dir.y < 0.) || !fRequireDownwards;

  bool end_is_fid = false;
  for (const geo::BoxBoundedGeo &g: fFiducialVolumes) {
    geo::Point_t end {t.end.x, t.end.y, t.end.z};

    bool is_contained = fCheckFiducialX ? g.ContainsPosition(end) : g.ContainsYZ(end.Y(), end.Z());

    if (is_contained) {
      end_is_fid = true;
      break;
    }
  }

  // Collection plane times need to be fiducial
  bool time_is_fid = \
    (t.hit_min_time_p2_tpcE < 0. || t.hit_min_time_p2_tpcE > fFidTickMin) &&
    (t.hit_max_time_p2_tpcE < 0. || t.hit_max_time_p2_tpcE < fFidTickMax) &&
    (t.hit_min_time_p2_tpcW < 0. || t.hit_min_time_p2_tpcW > fFidTickMin) &&
    (t.hit_max_time_p2_tpcW < 0. || t.hit_max_time_p2_tpcW < fFidTickMax);

  // compute the median dqdx of the last few cm -- using fMediandQdxRRMax
  std::vector<double> endp_dqdx;
  for (const sbn::TrackHitInfo &h: t.hits2) {
    if (h.oncalo && h.rr < fMediandQdxRRMax) endp_dqdx.push_back(h.dqdx);
  }
  double med_dqdx = -1;
  if (endp_dqdx.size()) {
    unsigned middle = endp_dqdx.size() / 2;
    std::nth_element(endp_dqdx.begin(), endp_dqdx.begin() + middle, endp_dqdx.end());
    med_dqdx = endp_dqdx[middle];

    // for even case
    if (endp_dqdx.size() % 2 == 0) {
      unsigned other_middle = middle - 1;
      std::nth_element(endp_dqdx.begin(), endp_dqdx.begin() + other_middle, endp_dqdx.end());
      med_dqdx = (med_dqdx + endp_dqdx[other_middle]) / 2.;
    }
  }

  bool valid_med_dqdx = ((med_dqdx > 0.) && (med_dqdx > fEndMediandQdxCut)) || (fEndMediandQdxCut < 0.);

  return downwards && end_is_fid && time_is_fid && valid_med_dqdx;
}

DEFINE_ART_CLASS_TOOL(TrackCaloSkimmerSelectStoppingTrack)

} // end namespace sbn
