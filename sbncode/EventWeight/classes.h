#include <map>
#include <string>
#include <vector>
#include "canvas/Persistency/Common/Wrapper.h"
#include "sbncode/EventWeight/MCEventWeight.h"

template class art::Wrapper<sbncode::evwgh::MCEventWeight>;
template class art::Wrapper<std::vector<sbncode::evwgh::MCEventWeight> >;
template class std::map<std::string, double>;
template class std::map<std::string, std::vector<double> >;

