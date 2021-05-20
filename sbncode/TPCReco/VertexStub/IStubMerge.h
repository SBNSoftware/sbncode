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
 *  @brief  IStubMerge interface class definiton
 */
class IStubMerge
{
public:
    /**
     *  @brief  Virtual Destructor
     */
    virtual ~IStubMerge() noexcept = default;

    virtual std::vector<sbn::StubInfo> Merge(const std::vector<sbn::StubInfo> &stubs,
       const geo::GeometryCore *geo, 
       const spacecharge::SpaceCharge *sce, 
       const detinfo::DetectorClocksData &dclock,
       const detinfo::DetectorPropertiesData &dprop) = 0;
};

} // namespace sbn
#endif

