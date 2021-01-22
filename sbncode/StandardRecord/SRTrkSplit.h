////////////////////////////////////////////////////////////////////////
// \file    SRTrkSplit.h
////////////////////////////////////////////////////////////////////////
#ifndef SRTRKSPLIT_H
#define SRTRKSPLIT_H

#include "SRVector3D.h"

namespace caf
{
  /// Representation of the reco momentum and PID a rb::Track from range
  class SRTrkSplit
    {
    public:

      SRTrkSplit();
      virtual ~SRTrkSplit();

      bool split;              ///< Whether this split exists
      SRVector3D locAtSplit;   ///< Location of the track at the split point
      SRVector3D dirAtSplit;   ///< Direction of the track at the split point
      int branch;              ///< The ID of the branch in the split
      int trunk;               ///< The ID of the trunk in the split
      float b_overlap;         ///< The overlap of the branch over the trunk in the split

      void setDefault();
    };

} // end namespace

#endif // SRTRKSPLIT_H
//////////////////////////////////////////////////////////////////////////////
