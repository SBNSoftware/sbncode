#ifndef _sbnumurecodata_CRTMatch_hh
#define _sbnumurecodata_CRTMatch_hh

#include "sbndcode/CRT/CRTProducts/CRTTrack.hh"
#include "sbndcode/CRT/CRTProducts/CRTHit.hh"

namespace numu {

struct CRTMatch {
  sbnd::crt::CRTTrack track;
  bool has_track_match;
  sbnd::crt::CRTHit hit;
  bool has_hit_match;
  double hit_distance;
  double match_time;
};
}

#endif
