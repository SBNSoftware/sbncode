#ifndef _sbncode_DynamicSelector_hh
#define _sbncode_DynamicSelector_hh

#include <functional>
#include <string>

#include "../Data/RecoEvent.h"

namespace numu {

void PrintClasses();
typedef std::function<bool (const numu::RecoTrack &, const numu::RecoEvent &)> TrackSelector;
TrackSelector Compile(const std::string &cutstr);

enum ROOTValueClass {
  RVTrack, RVTrueTrack, RVEvent
};

enum ROOTValueType {
  RVBool, RVFloat, RVInteger, RVUnsigned,
};

struct ROOTValue {
  ROOTValueClass base_class;
  int offset;
  int size;
  ROOTValueType type;
};

struct LiteralValue {
  bool is_float;
  bool is_valid;
  int data_int;
  float data_num;
};


ROOTValue MakeROOTValue(const std::string &name);
LiteralValue MakeROOTLiteral(const ROOTValue &value, const numu::RecoTrack &track, const numu::RecoEvent &event);

std::vector<std::string> MultiplyNames(const std::vector<std::vector<std::string>> &strings);
std::vector<numu::TrackSelector> MultiplySelectors(const std::vector<std::vector<std::string>> &track_selector_strings);

} // end namespace numu



#endif
