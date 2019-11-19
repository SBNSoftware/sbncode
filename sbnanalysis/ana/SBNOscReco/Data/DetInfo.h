#ifndef _sbnumurecodata_DetInfo_hh
#define _sbnumurecodata_DetInfo_hh

#include <vector>

#include "TVector3.h"

namespace numu {
/**
 * Holder for CRT hits 
 */
struct CRTHit {
  float time; //!< CRT Hit time
  bool has_coincidence; //!< Whether the hit requires a coincidence on two planes 
  float pes; //!< Number of PE's in hit
  TVector3 location; //!< Location of the hit
  TVector3 uncertainty; //!< Uncertainty on the location of the hit
};

struct FlashTriggerPrimitive {
  unsigned channel;
  std::vector<std::pair<int, unsigned>> triggers;
};

}
#endif
