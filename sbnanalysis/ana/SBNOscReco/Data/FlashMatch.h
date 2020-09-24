#ifndef _sbnumurecodata_FlashMatch_hh
#define _sbnumurecodata_FlashMatch_hh
namespace numu {

/**
 * Object for flash match to TPC track
 */
struct FlashMatch {
  bool present;
  float time; //!< Time of flash [us]
  float score; //!< score of flash matching
  unsigned pe; //!< number of photo electons in flash
};
}
#endif
