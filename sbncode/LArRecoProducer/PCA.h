#ifndef _sbncode_larrecoproducer_pca_h_
#define _sbncode_larrecoproducer_pca_h_ 

#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/PFParticle.h"

namespace sbnpca {

  float VecAngle(std::array<float, 2> A, std::array<float, 2> B);

  std::array<float, 2> HitVector(const recob::Hit &A, const geo::GeometryCore *geo, const detinfo::DetectorPropertiesData &dprop);
  
  std::array<float, 2> HitVector(const recob::Hit &A, const recob::Hit &B, const geo::GeometryCore *geo, const detinfo::DetectorPropertiesData &dprop);

  float HitDistance(const recob::Hit &A, const recob::Hit &B, const geo::GeometryCore *geo, const detinfo::DetectorPropertiesData &dprop);

  std::tuple<std::vector<art::Ptr<recob::Hit>>, std::vector<art::Ptr<recob::Hit>>, bool> GetNearestHits(
      const std::vector<art::Ptr<recob::Hit>> &hits, int ihit, float distance,
      const geo::GeometryCore *geo,
      const detinfo::DetectorPropertiesData &dprop);

  std::array<float, 2> HitPCAVec(const std::vector<art::Ptr<recob::Hit>> &hits, const recob::Hit &center,
               const geo::GeometryCore *geo, const detinfo::DetectorPropertiesData &dprop);

  std::array<float, 2> HitPCAVec(const std::vector<art::Ptr<recob::Hit>> &hits, const recob::Vertex &center,
               const geo::GeometryCore *geo, const detinfo::DetectorPropertiesData &dprop);

  std::array<float, 2> HitPCAVec(const std::vector<art::Ptr<recob::Hit>> &hits, std::array<float, 2> center,
               const geo::GeometryCore *geo, const detinfo::DetectorPropertiesData &dprop);

  std::array<float, 2> HitPCAEigen(const std::vector<art::Ptr<recob::Hit>> &hits, const art::Ptr<recob::Hit> &center,
      const geo::GeometryCore *geo, const detinfo::DetectorPropertiesData &dprop);

  float Vert2HitDistance(const recob::Hit &hit, const recob::Vertex &vert, const geo::GeometryCore *geo, const detinfo::DetectorPropertiesData &dprop);
  std::array<float, 2> Vert2HitCoord(const recob::Vertex &vert, const geo::PlaneID &planeID, const geo::GeometryCore *geo, const detinfo::DetectorPropertiesData &dprop);
}

#endif
