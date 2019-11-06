#ifndef _sbnumurecodata_DetInfo_hh
#define _sbnumurecodata_DetInfo_hh

namespace numu {
/**
 * Holder for CRT hits 
 */
struct CRTHit {
  float time; //!< CRT Hit time
  bool has_coincidence; //!< Whether the hit requires a coincidence on two planes 
  float pes; //!< Number of PE's in hit
};

}
#endif
