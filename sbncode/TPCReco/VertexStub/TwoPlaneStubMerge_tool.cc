/**
 *
 */

// Framework Includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Persistency/Common/PtrMaker.h"
#include "art/Utilities/ToolMacros.h"
#include "cetlib/cpu_timer.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// local includes
#include "sbncode/TPCReco/VertexStub/IStubMerge.h"
#include "sbncode/TPCReco/VertexStub/StubMergeAlgorithms.h"

// LArSoft includes
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcore/CoreUtils/ServiceUtil.h"

// std includes
#include <string>
#include <iostream>
#include <memory>

namespace sbn{
/**
 *  @brief  TwoPlaneStubMerge class definiton
 */
class TwoPlaneStubMerge : public sbn::IStubMerge {
public:
    /**
     *  @brief  Constructor
     */
    TwoPlaneStubMerge(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~TwoPlaneStubMerge();

    std::vector<sbn::StubInfo> Merge(const std::vector<sbn::StubInfo> &stubs,
       const geo::GeometryCore *geo, 
       const spacecharge::SpaceCharge *sce, 
       const detinfo::DetectorClocksData &dclock,
       const detinfo::DetectorPropertiesData &dprop) override;

    sbn::StubInfo MergeStubs(const sbn::StubInfo &A, const sbn::StubInfo &B,
      const geo::GeometryCore *geo,
      const spacecharge::SpaceCharge *sce,
      const detinfo::DetectorPropertiesData &dprop);

    bool SortStubs(const sbn::StubInfo &A, const sbn::StubInfo &B);

private:
  sbn::PlaneTransform fPlaneTransform;
  double fMaxMergeTOff;
  double fMaxMergeQOff;
  bool fRemoveDuplicateMerges;
  bool fSaveOldStubs;
};

TwoPlaneStubMerge::TwoPlaneStubMerge(fhicl::ParameterSet const &pset):
  fPlaneTransform(pset.get<fhicl::ParameterSet>("PlaneTransform")),
  fMaxMergeTOff(pset.get<double>("MaxMergeTOff")),
  fMaxMergeQOff(pset.get<double>("MaxMergeQOff")),
  fRemoveDuplicateMerges(pset.get<bool>("RemoveDuplicateMerges")),
  fSaveOldStubs(pset.get<bool>("SaveOldStubs"))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

TwoPlaneStubMerge::~TwoPlaneStubMerge()
{
}
    
// Return true if stub A is "better" than stub B
bool TwoPlaneStubMerge::SortStubs(const sbn::StubInfo &A, const sbn::StubInfo &B) {
  // Count the hits on the main stub
  int A_nhit = A.stub.CoreNHit();
  int B_nhit = B.stub.CoreNHit();

  // First -- try more hits
  if (A_nhit != B_nhit) return A_nhit > B_nhit;

  // Sort by desirability of plane: 2 better than 0 better than 1
  return ((A.vhit_hit->WireID().Plane + 1) % 3) < ((B.vhit_hit->WireID().Plane + 1) % 3);
}

sbn::StubInfo TwoPlaneStubMerge::MergeStubs(const sbn::StubInfo &A, const sbn::StubInfo &B,
    const geo::GeometryCore *geo,
    const spacecharge::SpaceCharge *sce,
    const detinfo::DetectorPropertiesData &dprop) {
  sbn::StubInfo ret; 

  // Metadata

  // Combine the hits
  ret.hits.insert(ret.hits.end(), A.hits.begin(), A.hits.end());
  ret.hits.insert(ret.hits.end(), A.hits.begin(), A.hits.end());

  // Figure out which is the "better" stub-plane -- probably the one with the lower pitch
  const sbn::StubInfo &best = SortStubs(A, B) ? A : B;
  const sbn::StubInfo &othr = SortStubs(A, B) ? B : A;

  ret.pfp = best.pfp;
  ret.vhit = best.vhit;
  ret.vhit_hit = best.vhit_hit;

  // The stub info
  // Vertex should be the same between the two
  ret.stub.vtx = A.stub.vtx;
  // The real thing -- combine the two positions to get a new endpoint
  geo::Point_t end = sbn::TwoStubEndPosition(fPlaneTransform, best, othr, geo, sce, dprop);
  ret.stub.end = TVector3(end.X(), end.Y(), end.Z());

  // In all the plane info stuff put the best stub first
  ret.stub.plane.push_back(best.stub.plane.front());
  ret.stub.plane.push_back(othr.stub.plane.front());

  ret.stub.hits.push_back(best.stub.hits.front());
  ret.stub.hits.push_back(othr.stub.hits.front());

  ret.stub.pitch.push_back(best.stub.pitch.front());
  ret.stub.pitch.push_back(othr.stub.pitch.front());

  ret.stub.trkpitch.push_back(best.stub.trkpitch.front());
  ret.stub.trkpitch.push_back(othr.stub.trkpitch.front());

  ret.stub.vtx_w.push_back(best.stub.vtx_w.front());
  ret.stub.vtx_w.push_back(othr.stub.vtx_w.front());

  ret.stub.hit_w.push_back(best.stub.hit_w.front());
  ret.stub.hit_w.push_back(othr.stub.hit_w.front());

  return ret;
}

std::vector<sbn::StubInfo> TwoPlaneStubMerge::Merge(const std::vector<sbn::StubInfo> &stubs,
    const geo::GeometryCore *geo, 
    const spacecharge::SpaceCharge *sce,
    const detinfo::DetectorClocksData &dclock,
    const detinfo::DetectorPropertiesData &dprop) {

  // Internal struct
  struct MergeInfo {
    float toff;
    float qoff;
    unsigned i;
    unsigned j;
  };

  std::vector<MergeInfo> merges;

  // Compute all of the possible cross-plane matches
  for (unsigned i_stub = 0; i_stub < stubs.size(); i_stub++) {
    for (unsigned j_stub = 0; j_stub < stubs.size(); j_stub++) {
      const sbn::StubInfo &A = stubs[i_stub];     
      const sbn::StubInfo &B = stubs[j_stub];     
      if (A.stub.plane.size() != 1 || B.stub.plane.size() != 1) continue; // Only match 1-plane-stubs

      if (A.stub.plane.front().Plane == B.stub.plane.front().Plane) continue; // Match across planes
      if (A.stub.plane.front().TPC != B.stub.plane.front().TPC) continue; // But should be in the same TPC

      float toff = sbn::StubTimeOffset(A, B, dclock, dprop);
      float qoff = sbn::StubChargeOffset(A, B);

      merges.push_back({toff, qoff, i_stub, j_stub});
    }
  }

  // Greedily do the best matches
  //
  // Sort by the time-offset -- this one is more precise
  std::sort(merges.begin(), merges.end(), [](auto const &lhs, auto const &rhs) {return lhs.toff < rhs.toff;});

  // keep track of which index stubs have been merged
  std::set<unsigned> hasmerged;

  // keep track of which index stubs should not be saved
  std::set<unsigned> dontsave;

  // Output list
  std::vector<sbn::StubInfo> ret;

  for (unsigned i_mrg = 0; i_mrg < merges.size(); i_mrg++) {
    const MergeInfo &thismrg = merges[i_mrg];
    if (thismrg.toff < fMaxMergeTOff && thismrg.qoff < fMaxMergeQOff) {
      if (!hasmerged.count(thismrg.i) && !hasmerged.count(thismrg.j)) { 
        ret.push_back(MergeStubs(stubs[thismrg.i], stubs[thismrg.j], geo, sce, dprop));
        hasmerged.insert(thismrg.i);
        hasmerged.insert(thismrg.j);
      }
      else if (!hasmerged.count(thismrg.i) && fRemoveDuplicateMerges) {
        dontsave.insert(thismrg.i);
      }
      else if (!hasmerged.count(thismrg.j) && fRemoveDuplicateMerges) {
        dontsave.insert(thismrg.j);
      }
    }
  }

  // Save the stubs which were not merged across planes and not marked to be removed
  for (unsigned i_stub = 0; i_stub < stubs.size(); i_stub++) {
    if (fSaveOldStubs || (!hasmerged.count(i_stub) && !dontsave.count(i_stub))) {
      ret.push_back(stubs[i_stub]);
    }
  }

  return ret;
}

DEFINE_ART_CLASS_TOOL(TwoPlaneStubMerge)

} // namespace sbn
