#include "GeoUtil.h"

std::vector<geo::BoxBoundedGeo> SBNRecoUtils::ActiveVolumes(const geo::GeometryCore *geometry) {
  std::vector<std::vector<geo::BoxBoundedGeo>> tpc_volumes = TPCVolumes(geometry);;
  std::vector<geo::BoxBoundedGeo> active_volumes;

  // make each cryostat volume a box enclosing all tpc volumes
  for (const std::vector<geo::BoxBoundedGeo> &tpcs: tpc_volumes) {
    double XMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinX() < rhs.MinX(); })->MinX();
    double YMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinY() < rhs.MinY(); })->MinY();
    double ZMin = std::min_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MinZ() < rhs.MinZ(); })->MinZ();
    
    double XMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxX() < rhs.MaxX(); })->MaxX();
    double YMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxY() < rhs.MaxY(); })->MaxY();
    double ZMax = std::max_element(tpcs.begin(), tpcs.end(), [](auto &lhs, auto &rhs) { return lhs.MaxZ() < rhs.MaxZ(); })->MaxZ();
    
    active_volumes.emplace_back(XMin, XMax, YMin, YMax, ZMin, ZMax);
  } 
  return active_volumes;
}

std::vector<std::vector<geo::BoxBoundedGeo>> SBNRecoUtils::TPCVolumes(const geo::GeometryCore *geometry) {
  std::vector<std::vector<geo::BoxBoundedGeo>> tpc_volumes;
  for (auto const &cryo: geometry->IterateCryostats()) {
    geo::GeometryCore::TPC_iterator iTPC = geometry->begin_TPC(cryo.ID()),
                                    tend = geometry->end_TPC(cryo.ID());
    std::vector<geo::BoxBoundedGeo> this_tpc_volumes;
    while (iTPC != tend) {
      geo::TPCGeo const& TPC = *iTPC;
      this_tpc_volumes.push_back(TPC.ActiveBoundingBox());
      iTPC++;
    }
     tpc_volumes.push_back(std::move(this_tpc_volumes));
  }
  return tpc_volumes;
}

