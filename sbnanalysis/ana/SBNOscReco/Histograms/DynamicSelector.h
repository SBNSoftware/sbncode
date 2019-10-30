#ifndef _sbncode_DynamicSelector_hh
#define _sbncode_DynamicSelector_hh

#include <functional>
#include <string>

#include "../Data/RecoEvent.h"

namespace numu {
  void PrintClasses();
  typedef std::function<bool (const numu::RecoTrack &, const numu::RecoEvent &)> TrackSelector;
  TrackSelector Compile(const std::string &cutstr);
}



#endif
