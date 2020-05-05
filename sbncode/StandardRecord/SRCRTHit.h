////////////////////////////////////////////////////////////////////////
// \file    SRCRTHit.h
////////////////////////////////////////////////////////////////////////
#ifndef SRCRTHIT_H
#define SRCRTHIT_H

#include "SRVector3D.h"

namespace caf
{
  /// Representation of the reco momentum and PID a recob::Track for 
  /// muon, pion, kaon, and proton assumptions 
  class SRCRTHit
    {
    public:

      SRCRTHit();
      virtual ~SRCRTHit() {}
      SRVector3D position;  // Position of CRT hit in detector coordinates [cm]
      SRVector3D position_err; // Error in position of CRT hit [cm]
      float time; // Time of CRT hit [us]
      void setDefault();
    };

} // end namespace

#endif // SRCRTHIT_H
//////////////////////////////////////////////////////////////////////////////
