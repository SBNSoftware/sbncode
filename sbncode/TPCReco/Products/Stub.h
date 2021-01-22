#ifndef sbncode_Stub_HH
#define sbncode_Stub_HH

#include <vector>
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "TVector3.h"

namespace sbn {
  class Stub {
  public:
    TVector3 vtx; //!< Interaction Vertex / Start of Stub [cm]
    TVector3 end; //!< End of Stub [cm]
    float charge; //!< Total charge of stub, corrected for elec. lifetime [#elec]
    float pitch; //!< Pitch of stub on each wire [cm]
    geo::PlaneID plane; //!< Plane of stub
    int nwire; //!< Number of wires/hits contributing charge to stub

    float wire_dist; //!< not used
  };
} // end namespace sbn

#endif
