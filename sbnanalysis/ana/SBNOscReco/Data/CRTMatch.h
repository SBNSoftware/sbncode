#ifndef _sbnumurecodata_CRTMatch_hh
#define _sbnumurecodata_CRTMatch_hh

#include "ana/SBNOscReco/Data/DetInfo.h"

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
  struct HitMatch {
    bool present; //!< Whether this CRTMatch has a matching hit
    float distance; //!< //!< Distance from projected track to CRT Hit. Nonsense if present is false.
    float time;
  };

  Track track; //!< CRT Track match
  HitMatch hit_match; //!< CRT Hit match
  CRTHit hit; // !< CRT hit object corresponsing to this one

};
}

#endif
