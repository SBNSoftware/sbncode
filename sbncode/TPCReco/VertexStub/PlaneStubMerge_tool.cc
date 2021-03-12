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
 *  @brief  PlaneStubMerge class definiton
 */
class PlaneStubMerge : public sbn::IStubMerge {
public:
    /**
     *  @brief  Constructor
     */
    PlaneStubMerge(fhicl::ParameterSet const &pset);

    /**
     *  @brief  Destructor
     */
    ~PlaneStubMerge();

    void Merge(std::vector<sbn::StubInfo> &stubs, const recob::Vertex &vertex, 
       const geo::GeometryCore *geo, 
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

//------------------------------------------------------------------------------------------------------------------------------------------

PlaneStubMerge::~PlaneStubMerge()
{
}

void PlaneStubMerge::Merge(std::vector<sbn::StubInfo> &stubs, const recob::Vertex &vertex, 
    const geo::GeometryCore *geo, 
    const spacecharge::SpaceCharge *sce,
    const detinfo::DetectorClocksData &dclock,
    const detinfo::DetectorPropertiesData &dprop) {
  for (unsigned i_stub = 0; i_stub < stubs.size(); i_stub++) {
    for (unsigned j_stub = 0; j_stub < stubs.size(); j_stub++) {
      if (i_stub == j_stub) continue;
      if (sbn::StubContains(stubs[i_stub], stubs[j_stub])) {
        if (sbn::StubDirectionDot(stubs[i_stub], stubs[j_stub], geo, dprop) > fStubDotCut) {
          // delete the smaller stub
          stubs.erase(stubs.begin() + j_stub);
        }
      }
    }
  }
}

DEFINE_ART_CLASS_TOOL(PlaneStubMerge)

} // namespace sbn
