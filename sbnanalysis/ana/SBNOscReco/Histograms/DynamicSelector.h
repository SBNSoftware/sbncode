#ifndef _sbncode_DynamicSelector_hh
#define _sbncode_DynamicSelector_hh

#include <functional>
#include <string>

#include "../uScript/value.h"
#include "../Data/RecoEvent.h"

namespace numu {

using TrackFunction = std::function<uscript::Value (const numu::RecoTrack *, const numu::TrueParticle *, const unsigned *)>;
using TrackSelector = std::function<bool (const numu::RecoTrack &, const numu::TrueParticle &, const unsigned &)>;

std::vector<std::string> MultiplyNames(const std::vector<std::vector<std::string>> &strings);
std::vector<TrackSelector> MultiplyTrackSelectors(const std::vector<std::vector<std::string>> &track_function_strings); 

} // end namespace numu



#endif
