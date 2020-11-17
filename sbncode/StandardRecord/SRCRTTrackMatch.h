////////////////////////////////////////////////////////////////////////
// \file    SRCRTTrackMatch.h
////////////////////////////////////////////////////////////////////////
#ifndef SRCRTTRACKMATCH_H
#define SRCRTTRACKMATCH_H

#include <limits>

namespace caf
{
  /// Matching information between a TPC Track and a CRT Track
  class SRCRTTrackMatch
    {
    public:

      SRCRTTrackMatch();
      virtual ~SRCRTTrackMatch() {}
      float time;  ///< Time of the CRT Track [us]
      float angle; ///< Relative angle between the TPC track and CRT track [rad]
    };

} // end namespace

#endif // SRCRTTRACKMATCH_H
//////////////////////////////////////////////////////////////////////////////
