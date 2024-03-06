/**
 *  @file   IStubMerge.h
 *
 *  @author grayputnam@uchicago.edu
 *
 */
#ifndef IStubMerge_h
#define IStubMerge_h

// Framework Includes
#include "sbnobj/Common/Reco/Stub.h"

#include "larcorealg/Geometry/GeometryCore.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesServiceStandard.h"
#include "lardataalg/DetectorInfo/DetectorClocksData.h"
#include "larevt/SpaceCharge/SpaceCharge.h"
#include "sbncode/TPCReco/VertexStub/StubMergeAlgorithms.h"

// cpp includes
#include <vector>

//------------------------------------------------------------------------------------------------------------------------------------------

namespace sbn {

/**
 *  @brief Abstract interface intended for art tools which take a list
 *         of stubs and return a new list with some of them merged.
 */
class IStubMerge
{
public:
    /**
     *  @brief  Virtual Destructor
     */
    virtual ~IStubMerge() noexcept = default;

    virtual std::vector<sbn::StubInfo> Merge(const std::vector<sbn::StubInfo> &stubs,
       const geo::WireReadoutGeom &channelMap,
       const spacecharge::SpaceCharge *sce, 
       const detinfo::DetectorClocksData &dclock,
       const detinfo::DetectorPropertiesData &dprop) = 0;
};

} // namespace sbn
#endif
