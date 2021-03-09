#ifndef StubMergeAlgorithms_HH
#define StubMergeAlgorithms_HH

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "sbncode/TPCReco/VertexStub/IStubMerge.h"

#include <memory>

namespace sbn {

bool StubContains(const sbn::StubInfo &A, const sbn::StubInfo &B);

float StubDirectionDot(const recob::Vertex &vertex, const sbn::StubInfo &A, const sbn::StubInfo &B, 
    const geo::GeometryCore *geo,
    const detinfo::DetectorPropertiesData &dprop);

} // end namespace sbn
#endif
