#ifndef _sbnumurecodata_CRTMatch_hh
#define _sbnumurecodata_CRTMatch_hh

#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"
#include "sbndcode/CRT/CRTProducts/CRTHit.hh"

namespace numu {

/**
* CRT matching object -- contains both a Hit and Track match
*/
struct CRTMatch {
  /**
  * Track Match
  */ 
  struct Track {
    bool present; //!< Whether this CRTMatch has a matching track
    float time; //!< Matching time [us] of track. T==0 is set to beam spill start time.
    float angle; //!< Angle between TPC track and CRT track
  };

  /**
  * Hit Match
  */
  struct Hit {
    bool present; //!< Whether this CRTMatch has a matching hit
    float time; //!< Matching time [us] of hit. T==0 is set to beam spill start time.
    float distance; //!< //!< Distance from projected track to CRT Hit. Nonsense if present is false.
  };

  Track track; //!< CRT Track match
  Hit hit; //!< CRT Hit match

};
}

#endif
