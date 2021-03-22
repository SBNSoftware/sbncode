#include <string>
#include <vector>
#include "sbnobj/Common/SBNEventWeight/EventWeightMap.h"
#include "WeightManager.h"

namespace sbn {
  namespace evwgh {

EventWeightMap WeightManager::Run(art::Event& e, const int inu) {
  EventWeightMap mcwgh;

  for (auto it=fWeightCalcMap.begin(); it!=fWeightCalcMap.end(); ++it) {
    const std::vector<float>& weights = it->second->GetWeight(e, inu);
    std::string wname = it->first + "_" + it->second->GetType();
    mcwgh.insert({ wname, weights });
  }

  return mcwgh;
}

  }  // namespace evwgh
}  // namespace sbn

