/**
 *  @file   IStubMerge.h
 *
 *  @author grayputnam@uchicago.edu
 *
 */
#ifndef IStubMerge_h
#define IStubMerge_h

// Framework Includes
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"

#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "sbnobj/Common/Reco/VertexHit.h"
#include "sbnobj/Common/Reco/Stub.h"

#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "lardataalg/DetectorInfo/DetectorPropertiesStandard.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataalg/DetectorInfo/DetectorClocksStandard.h"
#include "sbncode/TPCReco/VertexStub/StubMergeAlgorithms.h"

// cpp includes
#include <vector>
#include <array>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace sbn {

/**
 *  @brief  IStubMerge interface class definiton
 */
class IStubMerge
{
public:
    /**
     *  @brief  Virtual Destructor
     */
    virtual ~IStubMerge() noexcept = default;

    virtual void Merge(std::vector<sbn::StubInfo> &stubs, const recob::Vertex &vertex, 
       const geo::GeometryCore *geo, 
       const spacecharge::SpaceCharge *sce, 
       const detinfo::DetectorClocksData &dclock,
       const detinfo::DetectorPropertiesData &dprop) = 0;
};

} // namespace sbn
#endif

