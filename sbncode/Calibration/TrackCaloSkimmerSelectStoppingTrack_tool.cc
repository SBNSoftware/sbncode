// Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "fhiclcpp/ParameterSet.h"

#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larcore/Geometry/Geometry.h"
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

  double fExpFitRCut;
  double fFitResidualsCut;
  bool fRequireFit;
  double fEndMediandQdxCut;
  unsigned fNumberTimeSamples;

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

};

TrackCaloSkimmerSelectStoppingTrack::TrackCaloSkimmerSelectStoppingTrack(const fhicl::ParameterSet &p):
  fFVInsetMinX(p.get<double>("FVInsetMinX")),
  fFVInsetMaxX(p.get<double>("FVInsetMaxX")),
  fFVInsetMinY(p.get<double>("FVInsetMinY")),
  fFVInsetMaxY(p.get<double>("FVInsetMaxY")),
  fFVInsetMinZ(p.get<double>("FVInsetMinZ")),
  fFVInsetMaxZ(p.get<double>("FVInsetMaxZ")),
  fMinTimeTickInset(p.get<double>("MinTimeTickInset")),
  fMaxTimeTickInset(p.get<double>("MaxTimeTickInset")),
  fExpFitRCut(p.get<double>("ExpFitRCut")),
  fFitResidualsCut(p.get<double>("FitResidualsCut")),
  fRequireFit(p.get<bool>("RequireFit")),
  fEndMediandQdxCut(p.get<double>("EndMediandQdxCut")),
  fNumberTimeSamples(p.get<unsigned>("NumberTimeSamples"))
{
  // Get the fiducial volume info
  const geo::GeometryCore *geometry = lar::providerFrom<geo::Geometry>();

  // first the TPC volumes 
  for (auto const &cryo: geometry->IterateCryostats()) {
    geo::GeometryCore::TPC_iterator iTPC = geometry->begin_TPC(cryo.ID()),
                                    tend = geometry->end_TPC(cryo.ID());
    std::vector<geo::BoxBoundedGeo> this_tpc_volumes;
    while (iTPC != tend) {
      geo::TPCGeo const& TPC = *iTPC;
      this_tpc_volumes.push_back(TPC.ActiveBoundingBox());
      iTPC++;
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
  bool downwards = t.dir_y < 0.;

  bool end_is_fid = false;
  for (const geo::BoxBoundedGeo &g: fFiducialVolumes) {
    geo::Point_t end {t.end_x, t.end_y, t.end_z};

    if (g.ContainsPosition(end)) {
      end_is_fid = true;
      break;
    }
  }

  std::cout << "Track end: " << t.end_x << " " << t.end_y << " " << t.end_z << " is FID: " << end_is_fid << std::endl;

  // Collection plane times need to be fiducial
  bool time_is_fid = \
    (t.hit_min_time_p2_tpcE < 0. || t.hit_min_time_p2_tpcE > fFidTickMin) &&
    (t.hit_max_time_p2_tpcE < 0. || t.hit_max_time_p2_tpcE < fFidTickMax) &&
    (t.hit_min_time_p2_tpcW < 0. || t.hit_min_time_p2_tpcW > fFidTickMin) &&
    (t.hit_max_time_p2_tpcW < 0. || t.hit_max_time_p2_tpcW < fFidTickMax);

  std::cout << "Track times: " << t.hit_min_time_p2_tpcE << " " << t.hit_max_time_p2_tpcE << " " << t.hit_min_time_p2_tpcW << " " << t.hit_max_time_p2_tpcW << " is FID: " << time_is_fid << std::endl;

  // compute the median dqdx of the last 5 cm
  std::vector<double> endp_dqdx;
  for (const sbn::HitInfo &h: t.hits2) {
    if (h.oncalo && h.rr < 5.) endp_dqdx.push_back(h.dqdx);
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

  std::cout << "Track MED dqdx: " << med_dqdx << " is valid: " << valid_med_dqdx << std::endl;

  // Make sure the stopping fit worked
  bool valid_stopping_fit = t.n_fit_point > 2 || !fRequireFit;

  // Make a cut on the stopping fits
  bool valid_stopping_expR = (t.exp_fit_R < fExpFitRCut) || (fExpFitRCut < 0.) || !fRequireFit;

  double reduced_residual_offset = (t.const_fit_residuals - t.exp_fit_residuals) / t.n_fit_point;

  bool valid_stopping_residuals = (reduced_residual_offset < fFitResidualsCut) || (fFitResidualsCut < 0.) || !fRequireFit;

  std::cout << "Fit is valid: " << valid_stopping_fit << " " << valid_stopping_expR << " " << valid_stopping_residuals << std::endl;

  return downwards && end_is_fid && time_is_fid && valid_med_dqdx && valid_stopping_fit && valid_stopping_expR && valid_stopping_residuals;
}

DEFINE_ART_CLASS_TOOL(TrackCaloSkimmerSelectStoppingTrack)

} // end namespace sbn
