////////////////////////////////////////////////////////////////////////
// \file    SRCRTHitMatch.h
////////////////////////////////////////////////////////////////////////
#ifndef SRCRTHITMATCH_H
#define SRCRTHITMATCH_H

#include "SRCRTHit.h"
#include "SRVector3D.h"

namespace caf
{
  /// Representation of the reco momentum and PID a recob::Track for 
  /// muon, pion, kaon, and proton assumptions 
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
