#ifndef _sbnumurecodata_Mode_hh
#define _sbnumurecodata_Mode_hh

namespace numu {
/**
* Enum to hold each different typoe of reconstructed event
*/
enum InteractionMode {
  mCC = 0,  
  mNC = 1, 
  mCosmic = 2, 
  mOther = 3,
  mAll = 4
};

/**
 * Enum to hold different types of tracks.
 */
enum TrackMode {
  tmOther = -1,
  tmCosmic = 1,
  tmNeutrino = 2
};
}
#endif
