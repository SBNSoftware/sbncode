#include "PCA.h"

using namespace sbnpca;

float sbnpca::VecAngle(std::array<float, 2> A, std::array<float, 2> B) {
  float costh = (A[0] * B[0] + A[1] * B[1]) \
    / (sqrt(A[0]*A[0] + A[1] * A[1]) * sqrt(B[0]*B[0] + B[1] * B[1]));

  return acos(costh);
}

std::array<float, 2> sbnpca::HitVector(const recob::Hit &A, const geo::GeometryCore *geo, const detinfo::DetectorPropertiesData &dprop) {
  // get the wire distance between A and B
  float wire_distance = A.WireID().Wire;
  // convert to cm
  float wire_distance_cm = wire_distance * geo->WirePitch();

  // and the time difference
  float time_distance = A.PeakTime();
  // convert to cm
  float time_distance_cm = dprop.ConvertTicksToX(time_distance, A.WireID());

  return {wire_distance_cm, time_distance_cm};
}

std::array<float, 2> sbnpca::HitVector(const recob::Hit &A, const recob::Hit &B, const geo::GeometryCore *geo, const detinfo::DetectorPropertiesData &dprop) {
  // get each individual vec
  std::array<float, 2> vecA = HitVector(A, geo, dprop);
  std::array<float, 2> vecB = HitVector(B, geo, dprop);

  return {vecA[0] - vecB[0], vecA[1] - vecB[1]};
}

float sbnpca::HitDistance(const recob::Hit &A, const recob::Hit &B, const geo::GeometryCore *geo, const detinfo::DetectorPropertiesData &dprop) {
  std::array<float, 2> vec = HitVector(A, B, geo, dprop);
  return sqrt(vec[0]*vec[0] + vec[1]*vec[1]);
}

std::tuple<std::vector<art::Ptr<recob::Hit>>, std::vector<art::Ptr<recob::Hit>>, bool> sbnpca::GetNearestHits(
                                                 const std::vector<art::Ptr<recob::Hit>> &hits, int ihit, float distance,
                                                 const geo::GeometryCore *geo,
                                                 const detinfo::DetectorPropertiesData &dprop) {
  std::vector<art::Ptr<recob::Hit>> retlo;
  std::vector<art::Ptr<recob::Hit>> rethi;

  bool lo_complete = false;
  // pull in smaller ones
  for (int j = ihit-1; j >= 0; j--) {
    if (HitDistance(*hits[j], *hits[ihit], geo, dprop) > distance) {
      lo_complete = true;
      break;
    }
    retlo.push_back(hits[j]);
  }

  bool hi_complete = false;
  // pull in larger ones
  for (unsigned j = ihit+1; j < hits.size(); j++) {
    if (HitDistance(*hits[j], *hits[ihit], geo, dprop) > distance) {
      hi_complete = true;
      break;
    }
    rethi.push_back(hits[j]);
  }

  return {retlo, rethi, lo_complete && hi_complete};
}

std::array<float, 2> sbnpca::HitPCAVec(const std::vector<art::Ptr<recob::Hit>> &hits, const recob::Hit &center,
               const geo::GeometryCore *geo, const detinfo::DetectorPropertiesData &dprop) {

  return HitPCAVec(hits, HitVector(center, geo, dprop), geo, dprop);
}

std::array<float, 2> sbnpca::HitPCAVec(const std::vector<art::Ptr<recob::Hit>> &hits, const recob::Vertex &center,
               const geo::GeometryCore *geo, const detinfo::DetectorPropertiesData &dprop) {

  if (hits.size() == 0) return {-100., -100.};

  return HitPCAVec(hits, Vert2HitCoord(center, hits[0]->WireID(), geo, dprop), geo, dprop);
}

std::array<float, 2> sbnpca::HitPCAVec(const std::vector<art::Ptr<recob::Hit>> &hits, std::array<float, 2> center,
               const geo::GeometryCore *geo, const detinfo::DetectorPropertiesData &dprop) {

  if (hits.size() < 2) return {-100., -100.};

  std::array<float, 2> sum {};
  for (const art::Ptr<recob::Hit> &h: hits) {
    std::array<float, 2> vec = HitVector(*h, geo, dprop);
    sum[0] += vec[0] - center[0];
    sum[1] += vec[1] - center[1];
  }
  sum[0] = sum[0] / hits.size();
  sum[1] = sum[1] / hits.size();

  std::array<std::array<float, 2>, 2> scatter {};
  for (const art::Ptr<recob::Hit> &h: hits) {
    std::array<float, 2> vec = HitVector(*h, geo, dprop);
    vec[0] -= sum[0] - center[0];
    vec[1] -= sum[1] - center[1];

    scatter[0][0] += vec[0] * vec[0];
    scatter[0][1] += vec[0] * vec[1];
    scatter[1][0] += vec[1] * vec[0];
    scatter[1][1] += vec[1] * vec[1];
  }

  // first get the eigenvalues of the matrix
  float trace = scatter[0][0] + scatter[1][1];
  float det = scatter[0][0] * scatter[1][1] - scatter[0][1] * scatter[1][0];

  // this is always the max-eigenvalue
  float eigenP = (1. / 2.) * (trace + sqrt(trace*trace - 4 * det));
  // float eigenM = (1. / 2.) * (trace - sqrt(trace*trace - 4 * det));

  // and then the eigenvectors
  std::array<float, 2> ret {scatter[0][1], eigenP - scatter[0][0]};
  // std::array<float, 2> eigenVM {scatter[0][1], eigenM - scatter[0][0]};

  // make sure the sign is right
  if (sum[0] * ret[0] + sum[1] * ret[1] < 0.) {
    ret[0] = -ret[0];
    ret[1] = -ret[1];
  }

  return ret;

}

std::array<float, 2> sbnpca::HitPCAEigen(const std::vector<art::Ptr<recob::Hit>> &hits, const art::Ptr<recob::Hit> &center,
               const geo::GeometryCore *geo, const detinfo::DetectorPropertiesData &dprop) {

  std::array<float, 2> sum {};
  for (const art::Ptr<recob::Hit> &h: hits) {
    std::array<float, 2> vec = HitVector(*h, *center, geo, dprop);
    sum[0] += vec[0];
    sum[1] += vec[1];
  }
  sum[0] = sum[0] / hits.size();
  sum[1] = sum[1] / hits.size();

  std::array<std::array<float, 2>, 2> scatter {};
  for (const art::Ptr<recob::Hit> &h: hits) {
    std::array<float, 2> vec = HitVector(*h, *center, geo, dprop);
    vec[0] -= sum[0];
    vec[1] -= sum[1];

    scatter[0][0] += vec[0] * vec[0];
    scatter[0][1] += vec[0] * vec[1];
    scatter[1][0] += vec[1] * vec[0];
    scatter[1][1] += vec[1] * vec[1];
  }

  // first get the eigenvalues of the matrix
  float trace = scatter[0][0] + scatter[1][1];
  float det = scatter[0][0] * scatter[1][1] - scatter[0][1] * scatter[1][0];

  float eigenP = (1. / 2.) * (trace + sqrt(trace*trace - 4 * det));
  float eigenM = (1. / 2.) * (trace - sqrt(trace*trace - 4 * det));

  return {eigenP, eigenM};
}

std::array<float, 2> sbnpca::Vert2HitCoord(const recob::Vertex &vert, const geo::PlaneID &planeID, 
    const geo::GeometryCore *geo, const detinfo::DetectorPropertiesData &dprop) {

  float vert_wire_coord = geo->WireCoordinate(vert.position().Y(), vert.position().Z(), planeID) * geo->WirePitch();
  float vert_time_coord = vert.position().X();

  return {vert_wire_coord, vert_time_coord};
}

float sbnpca::Vert2HitDistance(const recob::Hit &hit, const recob::Vertex &vert, const geo::GeometryCore *geo, const detinfo::DetectorPropertiesData &dprop) {
  float vert_wire_coord = geo->WireCoordinate(vert.position().Y(), vert.position().Z(), hit.WireID()) * geo->WirePitch();
  float hit_wire_coord = hit.WireID().Wire * geo->WirePitch();

  float vert_time_coord = vert.position().X();
  float hit_time_coord = dprop.ConvertTicksToX(hit.PeakTime(), hit.WireID());

  return sqrt((vert_wire_coord - hit_wire_coord) * (vert_wire_coord - hit_wire_coord) +
    (vert_time_coord - hit_time_coord) * (vert_time_coord - hit_time_coord));
}


