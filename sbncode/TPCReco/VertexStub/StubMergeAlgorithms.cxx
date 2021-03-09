#include "StubMergeAlgorithms.h"

// Check if A contains B
bool sbn::StubContains(const sbn::StubInfo &A, const sbn::StubInfo &B) {
  if (A.vhit_hit->WireID().Plane != B.vhit_hit->WireID().Plane) return false;

  for (const art::Ptr<recob::Hit> &h: A.hits) {
    if (h.key() == B.vhit_hit.key()) return true;
  }
  return false;
}

float sbn::StubDirectionDot(const recob::Vertex &vertex, const sbn::StubInfo &A, const sbn::StubInfo &B, 
    const geo::GeometryCore *geo,
    const detinfo::DetectorPropertiesData &dprop) {
  
  // project the vertex onto the wireplane
  float vert_w = geo->WireCoordinate(vertex.position(), A.vhit_hit->WireID()) * geo->WirePitch();
  float vert_x = vertex.position().x();

  float Ahit_w = A.vhit_hit->WireID().Wire * geo->WirePitch() - vert_w;
  float Ahit_x = dprop.ConvertTicksToX(A.vhit_hit->PeakTime(), A.vhit_hit->WireID()) - vert_x;

  float Bhit_w = B.vhit_hit->WireID().Wire * geo->WirePitch() - vert_w;
  float Bhit_x = dprop.ConvertTicksToX(B.vhit_hit->PeakTime(), B.vhit_hit->WireID()) - vert_x;

  float Amag = sqrt(Ahit_w*Ahit_w + Ahit_x*Ahit_x);
  float Bmag = sqrt(Bhit_w*Bhit_w + Bhit_x*Bhit_x);

  return(Ahit_w*Bhit_w + Ahit_x*Bhit_x) / (Amag*Bmag);
}

