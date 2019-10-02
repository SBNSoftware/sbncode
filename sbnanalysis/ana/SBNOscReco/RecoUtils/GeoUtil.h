#ifndef _sbnanalysis_numu_GeoUtil_hh_
#define _sbnanalysis_numu_GeoUtil_hh_
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larcorealg/Geometry/GeometryCore.h"
namespace SBNRecoUtils {
  std::vector<geo::BoxBoundedGeo> ActiveVolumes(const geo::GeometryCore *geometry);
  std::vector<std::vector<geo::BoxBoundedGeo>> TPCVolumes(const geo::GeometryCore *geometry);
}
#endif
