////////////////////////////////////////////////////////////////////////
// \file    SRCRTTrack.h
////////////////////////////////////////////////////////////////////////
#ifndef SRCRTTRACK_H
#define SRCRTTRACK_H

#include "sbncode/StandardRecord/SRCRTHit.h"

namespace caf
{
  class SRCRTTrack
    {
    public:

      SRCRTTrack();
      virtual ~SRCRTTrack() {}

      SRCRTHit hita; //!< First hit in CRT track
      SRCRTHit hitb; //!< Second git in CRT Track
      float time;    //!< Combined time of CRT Track
    };

} // end namespace

#endif // SRCRTTRACK_H
//////////////////////////////////////////////////////////////////////////////
