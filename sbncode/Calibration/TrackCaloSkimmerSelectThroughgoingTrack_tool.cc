#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include "larcore/Geometry/Geometry.h"
#include "larcore/CoreUtils/ServiceUtil.h"
#include "larcorealg/Geometry/GeometryCore.h"

#include "ITCSSelectionTool.h"

namespace sbn {

  class TrackCaloSkimmerSelectThroughgoingTrack: public ITCSSelectionTool {
  public:

    TrackCaloSkimmerSelectThroughgoingTrack(const fhicl::ParameterSet &p);
    ~TrackCaloSkimmerSelectThroughgoingTrack() {}

    bool Select(const TrackInfo &t) override;

  private:
    double fFVInsetMinX;
    double fFVInsetMaxX;
    double fFVInsetMinY;
    double fFVInsetMaxY;
    double fFVInsetMinZ;
    double fFVInsetMaxZ;
    bool   fCheckFiducialX;

    std::vector<geo::BoxBoundedGeo> fTPCVolumes;
    geo::BoxBoundedGeo fActiveVolume, fFiducialVolume;
  };

  TrackCaloSkimmerSelectThroughgoingTrack::TrackCaloSkimmerSelectThroughgoingTrack(const fhicl::ParameterSet &p):
    ITCSSelectionTool(p),
    fFVInsetMinX(p.get<double>("FVInsetMinX")),
    fFVInsetMaxX(p.get<double>("FVInsetMaxX")),
    fFVInsetMinY(p.get<double>("FVInsetMinY")),
    fFVInsetMaxY(p.get<double>("FVInsetMaxY")),
    fFVInsetMinZ(p.get<double>("FVInsetMinZ")),
    fFVInsetMaxZ(p.get<double>("FVInsetMaxZ")),
    fCheckFiducialX(p.get<bool>("CheckFiducialX"))
  {
    const geo::GeometryCore *geometry = lar::providerFrom<geo::Geometry>();

    for(auto const &cryoid: geometry->Iterate<geo::CryostatID>())
      {
        for(auto const& TPC : geometry->Iterate<geo::TPCGeo>(cryoid))
          fTPCVolumes.push_back(TPC.ActiveBoundingBox());
      }

    double XMin = std::min_element(fTPCVolumes.begin(), fTPCVolumes.end(), [](auto &lhs, auto &rhs) { return lhs.MinX() < rhs.MinX(); })->MinX();
    double YMin = std::min_element(fTPCVolumes.begin(), fTPCVolumes.end(), [](auto &lhs, auto &rhs) { return lhs.MinY() < rhs.MinY(); })->MinY();
    double ZMin = std::min_element(fTPCVolumes.begin(), fTPCVolumes.end(), [](auto &lhs, auto &rhs) { return lhs.MinZ() < rhs.MinZ(); })->MinZ();

    double XMax = std::max_element(fTPCVolumes.begin(), fTPCVolumes.end(), [](auto &lhs, auto &rhs) { return lhs.MaxX() < rhs.MaxX(); })->MaxX();
    double YMax = std::max_element(fTPCVolumes.begin(), fTPCVolumes.end(), [](auto &lhs, auto &rhs) { return lhs.MaxY() < rhs.MaxY(); })->MaxY();
    double ZMax = std::max_element(fTPCVolumes.begin(), fTPCVolumes.end(), [](auto &lhs, auto &rhs) { return lhs.MaxZ() < rhs.MaxZ(); })->MaxZ();

    fActiveVolume = geo::BoxBoundedGeo(XMin, XMax, YMin, YMax, ZMin, ZMax);

    fFiducialVolume = geo::BoxBoundedGeo(fActiveVolume.MinX() + fFVInsetMinX, fActiveVolume.MaxX() - fFVInsetMaxX,
                                         fActiveVolume.MinY() + fFVInsetMinY, fActiveVolume.MaxY() - fFVInsetMaxY,
                                         fActiveVolume.MinZ() + fFVInsetMinZ, fActiveVolume.MaxZ() - fFVInsetMaxZ);
  }

  bool TrackCaloSkimmerSelectThroughgoingTrack::Select(const TrackInfo &t)
  {
    geo::Point_t start {t.start.x, t.start.y, t.start.z};
    bool start_is_non_fid = fCheckFiducialX ? !fFiducialVolume.ContainsPosition(start) : !fFiducialVolume.ContainsYZ(start.Y(), start.Z());

    geo::Point_t end {t.end.x, t.end.y, t.end.z};
    bool end_is_non_fid = fCheckFiducialX ? !fFiducialVolume.ContainsPosition(end) : !fFiducialVolume.ContainsYZ(end.Y(), end.Z());

    return start_is_non_fid && end_is_non_fid;
  }

  DEFINE_ART_CLASS_TOOL(TrackCaloSkimmerSelectThroughgoingTrack)

}
