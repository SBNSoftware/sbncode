////////////////////////////////////////////////////////////////////////
// \file    SRCRTHitMatch.h
////////////////////////////////////////////////////////////////////////
#ifndef SRCRTHITMATCH_H
#define SRCRTHITMATCH_H

#include "sbncode/StandardRecord/SRCRTHit.h"
#include "sbncode/StandardRecord/SRVector3D.h"

namespace caf
{
  /// Matching information between a TPC Track and a CRT Hit
  class SRCRTHitMatch
    {
    public:

      SRCRTHitMatch();
      virtual ~SRCRTHitMatch() {}
      SRCRTHit hit; ///< The CRT hit
      float distance; ///< Distance from the projected TPC track to the CRT hit [cm]
    };

} // end namespace

#endif // SRCRTHITMATCH_H
//////////////////////////////////////////////////////////////////////////////
