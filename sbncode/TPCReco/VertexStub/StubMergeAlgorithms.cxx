#include "StubMergeAlgorithms.h"

geo::Point_t sbn::GetLocation(const spacecharge::SpaceCharge *sce, geo::Point_t loc_w, geo::TPCID TPC, float xsign) {
  if (sce && sce->EnableCalSpatialSCE()) {
    geo::Vector_t offset = sce->GetCalPosOffsets(loc_w, TPC.TPC);
    loc_w.SetX(loc_w.X() + xsign * offset.X());
    loc_w.SetY(loc_w.Y() + offset.Y());
    loc_w.SetZ(loc_w.Z() + offset.Z());
  }

  return loc_w;
}

double sbn::GetEfield(const detinfo::DetectorPropertiesData &dprop, const spacecharge::SpaceCharge *sce, geo::Point_t loc, geo::TPCID TPC, bool correct_loc_sce, float xsign) {

  double EField = dprop.Efield();
  if (sce && sce->EnableSimEfieldSCE()) {
    if (correct_loc_sce) loc = sbn::GetLocation(sce, loc, TPC, xsign);

    // Gets relative E field Distortions
    geo::Vector_t EFieldOffsets = sce->GetEfieldOffsets(loc);
    // Add 1 in X direction as this is the direction of the drift field
    EFieldOffsets = EFieldOffsets + geo::Xaxis();
    // Convert to Absolute E Field from relative
    EFieldOffsets = EField * EFieldOffsets;
    // We only care about the magnitude for recombination
    EField = EFieldOffsets.r();
  }

  return EField;
}

geo::Point_t sbn::GetLocationAtWires(const spacecharge::SpaceCharge *sce, const geo::GeometryCore *geo, geo::Point_t loc, geo::TPCID TPC, float xsign) { 
  if (sce && sce->EnableCalSpatialSCE()) {
    // Returned X is the drift -- multiply by the drift direction to undo this
    int corr = geo->TPC(TPC).DriftDir()[0];

    // for some reason, one needs to flip the sign of the x-direction when correcting for field distortion
    geo::Vector_t offset = sce->GetPosOffsets(loc);
    loc.SetX(loc.X() + corr * xsign * offset.X());
    loc.SetY(loc.Y() + offset.Y());
    loc.SetZ(loc.Z() + offset.Z());
  }

  return loc;
}

double sbn::GetPitch(
    const geo::GeometryCore *geo, const spacecharge::SpaceCharge *sce, 
    geo::Point_t loc, geo::Vector_t dir, 
    geo::View_t view, geo::TPCID tpc, 
    bool correct_sce, bool track_is_sce_corrected, float xsign) {

  double angleToVert = geo->WireAngleToVertical(view, tpc) - 0.5*::util::pi<>();

  geo::Vector_t dir_w;

  // "dir_w" should be the direction that the wires see. If the track already has the field
  // distortion corrections applied, then we need to de-apply them to get the direction as
  // seen by the wire planes
  if (sce && sce->EnableCalSpatialSCE() && correct_sce && track_is_sce_corrected) {
    // compute the dir of the track trajectory
    geo::Point_t loc_mdx = loc - dir * (geo->WirePitch(view) / 2.);
    geo::Point_t loc_pdx = loc + dir * (geo->WirePitch(view) / 2.);

    // map to the wires
    loc_mdx = GetLocationAtWires(sce, geo, loc_mdx, tpc, xsign);
    loc_pdx = GetLocationAtWires(sce, geo, loc_pdx, tpc, xsign);

    dir_w = (loc_pdx - loc_mdx).Unit(); 
  }
  // If there is no space charge or the track is not yet corrected, then the dir
  // is the track is what we want
  else {
    dir_w = dir;
  }

  double cosgamma = std::abs(std::sin(angleToVert)*dir_w.Y() + std::cos(angleToVert)*dir_w.Z());
  double pitch;
  if (cosgamma) {
    pitch = geo->WirePitch(view)/cosgamma;
  }
  else {
    pitch = 0.;
  }

  // now take the pitch computed on the wires and correct it back to the particle trajectory
  geo::Point_t loc_w = loc;
  if (correct_sce && track_is_sce_corrected) {
    loc_w = sbn::GetLocationAtWires(sce, geo, loc, tpc, xsign);
  }

  geo::Point_t locw_traj = sbn::GetLocation(sce, loc_w, tpc, xsign);
  geo::Point_t locw_pdx_traj = sbn::GetLocation(sce, loc_w + pitch * dir_w, tpc, xsign);

  pitch = (locw_traj - locw_pdx_traj).R();

  return pitch;
}


// Check if A contains B
bool sbn::StubContains(const sbn::StubInfo &A, const sbn::StubInfo &B) {
  if (A.vhit_hit->WireID().Plane != B.vhit_hit->WireID().Plane) return false;

  for (const art::Ptr<recob::Hit> &h: A.hits) {
    if (A.stub.OnCore(h->WireID())) {
      if (h.key() == B.vhit_hit.key()) {
        return true;
      }
    }
  }
  return false;
}

float sbn::StubDirectionDot(const sbn::StubInfo &A, const sbn::StubInfo &B, 
    const geo::GeometryCore *geo,
    const detinfo::DetectorPropertiesData &dprop) {

  geo::Point_t vertex(A.stub.vtx);
  
  // project the vertex onto the wireplane
  float vert_w = geo->WireCoordinate(vertex, A.vhit_hit->WireID()) * geo->WirePitch();
  float vert_x = vertex.x();

  float Ahit_w = A.vhit_hit->WireID().Wire * geo->WirePitch() - vert_w;
  float Ahit_x = dprop.ConvertTicksToX(A.vhit_hit->PeakTime(), A.vhit_hit->WireID()) - vert_x;

  float Bhit_w = B.vhit_hit->WireID().Wire * geo->WirePitch() - vert_w;
  float Bhit_x = dprop.ConvertTicksToX(B.vhit_hit->PeakTime(), B.vhit_hit->WireID()) - vert_x;

  float Amag = sqrt(Ahit_w*Ahit_w + Ahit_x*Ahit_x);
  float Bmag = sqrt(Bhit_w*Bhit_w + Bhit_x*Bhit_x);

  return(Ahit_w*Bhit_w + Ahit_x*Bhit_x) / (Amag*Bmag);
}

float sbn::StubTimeOffset(const sbn::StubInfo &A, const sbn::StubInfo &B, 
    const detinfo::DetectorClocksData &dclock,
    const detinfo::DetectorPropertiesData &dprop) {
  // Convert the ticks to an X-point to undo tick-offsets between planes. And then convert to ticks.
  return abs(dprop.ConvertTicksToX(A.vhit_hit->PeakTime(), A.vhit_hit->WireID()) - dprop.ConvertTicksToX(B.vhit_hit->PeakTime(), B.vhit_hit->WireID())) / dprop.DriftVelocity() / dclock.TPCClock().TickPeriod();
}

geo::Point_t sbn::TwoStubEndPosition(const sbn::PlaneTransform &T, const sbn::StubInfo &A, const sbn::StubInfo &B,
    const geo::GeometryCore *geo,
    const spacecharge::SpaceCharge *sce,
    const detinfo::DetectorPropertiesData &dprop) {

  static constexpr bool VERBOSE = false;

  // look up the wire-coordinate of both wires
  double A_w = T.WireCoordinate(A.vhit_hit->WireID());
  double B_w = T.WireCoordinate(B.vhit_hit->WireID());

  // Use this to get y, z
  double y = T.TwoPlaneToY(A.vhit_hit->WireID(), A_w, B.vhit_hit->WireID(), B_w);
  double z = T.TwoPlaneToZ(A.vhit_hit->WireID(), A_w, B.vhit_hit->WireID(), B_w);

  // just average the x-pos
  double x = (dprop.ConvertTicksToX(A.vhit_hit->PeakTime(), A.vhit_hit->WireID()) + dprop.ConvertTicksToX(B.vhit_hit->PeakTime(), B.vhit_hit->WireID())) / 2.;

  geo::Point_t pos(x, y, z);

  // Correct for space charge
  pos = GetLocation(sce, pos, A.vhit_hit->WireID());

  if (VERBOSE) {
    std::cout << "A View: " << geo->View(A.vhit_hit->WireID()) << " Angle to Vertical: " << geo->WireAngleToVertical(geo->View(A.vhit_hit->WireID()), A.vhit_hit->WireID()) << std::endl;
    std::cout << "B View: " << geo->View(B.vhit_hit->WireID()) << " Angle to Vertical: " << geo->WireAngleToVertical(geo->View(B.vhit_hit->WireID()), B.vhit_hit->WireID()) << std::endl;
    std::cout << "A Wire Pos Y: " << geo->Wire(A.vhit_hit->WireID()).GetCenter().Y() << " Z: " << geo->Wire(A.vhit_hit->WireID()).GetCenter().Z() << std::endl;
    std::cout << "B Wire Pos Y: " << geo->Wire(B.vhit_hit->WireID()).GetCenter().Y() << " Z: " << geo->Wire(B.vhit_hit->WireID()).GetCenter().Z() << std::endl;
    
    std::cout << "A Wire: " << A.vhit_hit->WireID().Wire << " Plane: " << A.vhit_hit->WireID().Plane << " TPC: " << A.vhit_hit->WireID().TPC << " Coord: " << A_w << std::endl;
    std::cout << "B Wire: " << B.vhit_hit->WireID().Wire << " Plane: " << B.vhit_hit->WireID().Plane << " TPC: " << B.vhit_hit->WireID().TPC << " Coord: " << B_w << std::endl;
    
    std::cout << "Output Y: " << y << " Z: " << z << std::endl;
    
    std::cout << "A Wire coord: " << geo->WireCoordinate(y, z, A.vhit_hit->WireID()) << std::endl;
    std::cout << "B Wire coord: " << geo->WireCoordinate(y, z, B.vhit_hit->WireID()) << std::endl;
    
    std::cout << "A X: " << dprop.ConvertTicksToX(A.vhit_hit->PeakTime(), A.vhit_hit->WireID()) << std::endl;
    std::cout << "B X: " << dprop.ConvertTicksToX(B.vhit_hit->PeakTime(), B.vhit_hit->WireID()) << std::endl;
    
    std::cout << "Output X: " << x << std::endl;
    
    
    std::cout << "Final x: " << pos.x() << " y: " << pos.y() << " z: " << pos.z() << std::endl;
  }

  return pos;
}

float sbn::StubChargeOffset(const sbn::StubInfo &A, const sbn::StubInfo &B) {
  assert(A.stub.hits.size() == 1 && B.stub.hits.size() == 1);

  // Only count charge on the main stub
  float ACharge = A.stub.CoreCharge();
  float BCharge = B.stub.CoreCharge();

  return abs(ACharge - BCharge);
}

float sbn::StubPeakChargeOffset(const sbn::StubInfo &A, const sbn::StubInfo &B) {
  return abs(A.vhit->charge - B.vhit->charge);
}

float sbn::StubPeakdQdxOffset(const sbn::PlaneTransform &T, const sbn::StubInfo &A, const sbn::StubInfo &B,
    const geo::GeometryCore *geo,
    const spacecharge::SpaceCharge *sce,
    const detinfo::DetectorPropertiesData &dprop) {

  TVector3 vtx = A.stub.vtx;
  geo::Point_t end = sbn::TwoStubEndPosition(T, A, B, geo, sce, dprop); 
  TVector3 end_v(end.x(), end.y(), end.z());

  geo::Vector_t dir((end_v-vtx).Unit());

  float Apitch = GetPitch(geo, sce, end, dir, A.vhit_hit->View(), A.vhit_hit->WireID(), true, true);
  float Bpitch = GetPitch(geo, sce, end, dir, B.vhit_hit->View(), B.vhit_hit->WireID(), true, true);

  return abs(A.vhit->charge / Apitch - B.vhit->charge / Bpitch);
}

