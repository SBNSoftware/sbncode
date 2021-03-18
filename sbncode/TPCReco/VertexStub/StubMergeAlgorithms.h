#ifndef StubMergeAlgorithms_HH
#define StubMergeAlgorithms_HH

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorClocksStandard.h"

#include "larevt/SpaceCharge/SpaceCharge.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "sbnobj/Common/Reco/VertexHit.h"
#include "sbnobj/Common/Reco/Stub.h"

#include "sbncode/TPCReco/VertexStub/PlaneTransform.h"

#include <memory>

namespace sbn {

struct StubInfo {
  sbn::Stub stub;
  art::Ptr<recob::PFParticle> pfp;
  std::vector<art::Ptr<recob::Hit>> hits;
  art::Ptr<sbn::VertexHit> vhit;
  art::Ptr<recob::Hit> vhit_hit;
};

// Helpers
double GetPitch(
    const geo::GeometryCore *geo, const spacecharge::SpaceCharge *sce, 
    geo::Point_t loc, geo::Vector_t dir, 
    geo::View_t view, geo::TPCID tpc, 
    bool correct_sce, bool track_is_sce_corrected, float xsign=1.);

geo::Point_t GetLocation(const spacecharge::SpaceCharge *sce, geo::Point_t loc_w, unsigned TPC, float xsign=1.);
double GetEfield(const detinfo::DetectorPropertiesData &dprop, const spacecharge::SpaceCharge *sce, geo::Point_t loc, unsigned TPC, bool correct_loc_sce, float xsign=1.);
geo::Point_t GetLocationAtWires(const spacecharge::SpaceCharge *sce, geo::Point_t loc, float xsign=1.);

// For matching stubs on the same plane
bool StubContains(const sbn::StubInfo &A, const sbn::StubInfo &B);

float StubDirectionDot(const sbn::StubInfo &A, const sbn::StubInfo &B, 
    const geo::GeometryCore *geo,
    const detinfo::DetectorPropertiesData &dprop);

// For matching stubs across planes
float StubTimeOffset(const sbn::StubInfo &A, const sbn::StubInfo &B, 
    const detinfo::DetectorClocksData &dclock,
    const detinfo::DetectorPropertiesData &dprop);

geo::Point_t TwoStubEndPosition(const sbn::PlaneTransform &T, const sbn::StubInfo &A, const sbn::StubInfo &B,
    const geo::GeometryCore *geo,
    const spacecharge::SpaceCharge *sce,
    const detinfo::DetectorPropertiesData &dprop);

float StubChargeOffset(const sbn::StubInfo &A, const sbn::StubInfo &B);
float StubPeakChargeOffset(const sbn::StubInfo &A, const sbn::StubInfo &B);
float StubPeakdQdxOffset(const sbn::PlaneTransform &T, const sbn::StubInfo &A, const sbn::StubInfo &B,
    const geo::GeometryCore *geo,
    const spacecharge::SpaceCharge *sce,
    const detinfo::DetectorPropertiesData &dprop);

} // end namespace sbn
#endif
