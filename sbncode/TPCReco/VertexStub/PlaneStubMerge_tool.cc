/**
 *
 */

// Framework Includes
#include "art/Utilities/ToolMacros.h"
#include "cetlib/cpu_timer.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// local includes
#include "sbncode/TPCReco/VertexStub/IStubMerge.h"
#include "sbncode/TPCReco/VertexStub/StubMergeAlgorithms.h"

// LArSoft includes
#include "larcorealg/Geometry/GeometryCore.h"

// std includes
#include <memory>

namespace sbn{
/**
 *  @brief  Art tool for merging stubs on the same plane.
 *
 * Always merges into the larger stub.
 */
class PlaneStubMerge : public sbn::IStubMerge {
public:
    /**
     *  @brief  Constructor
     */
    PlaneStubMerge(fhicl::ParameterSet const &pset);

    std::vector<sbn::StubInfo> Merge(const std::vector<sbn::StubInfo> &stubs,
       const geo::WireReadoutGeom &wireReadout,
       const spacecharge::SpaceCharge *sce, 
       const detinfo::DetectorClocksData &dclock,
       const detinfo::DetectorPropertiesData &dprop) override;

private:
  double fStubDotCut;
};

PlaneStubMerge::PlaneStubMerge(fhicl::ParameterSet const &pset):
  fStubDotCut(pset.get<double>("StubDotCut"))
{
}

std::vector<sbn::StubInfo> PlaneStubMerge::Merge(const std::vector<sbn::StubInfo> &stubs,
    const geo::WireReadoutGeom &wireReadout,
    const spacecharge::SpaceCharge *sce,
    const detinfo::DetectorClocksData &dclock,
    const detinfo::DetectorPropertiesData &dprop) {

  std::vector<bool> toerase(stubs.size(), false);
  for (unsigned i_stub = 0; i_stub < stubs.size(); i_stub++) {
    for (unsigned j_stub = 0; j_stub < stubs.size(); j_stub++) {
      if (i_stub == j_stub) continue;
      if (toerase[j_stub]) continue; // Already being erased
      if (sbn::StubContains(stubs[i_stub], stubs[j_stub])) {
        if (sbn::StubDirectionDot(stubs[i_stub], stubs[j_stub], wireReadout, dprop) > fStubDotCut) {
          toerase[j_stub] = true;
        }
      }
    }
  }

  std::vector<sbn::StubInfo> ret;
  // Keep only good stubs
  for (unsigned i = 0; i < stubs.size(); i++) {
    if (!toerase[i]) ret.push_back(stubs[i]);
  }

  return ret;
}

DEFINE_ART_CLASS_TOOL(PlaneStubMerge)

} // namespace sbn
