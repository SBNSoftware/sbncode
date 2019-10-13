#ifndef _sbnumurecodata_FlashMatch_hh
#define _sbnumurecodata_FlashMatch_hh
namespace numu {

/**
 * Object for flash match to TPC track
 */
struct FlashMatch {
  float match_time; //!< [us] PE weighted average of all matching hit times.
  float match_time_first; //!< [us] First matching hit time.
  float match_time_width; //!< [us] Width of hit match. Currently not set.
};
}
#endif
