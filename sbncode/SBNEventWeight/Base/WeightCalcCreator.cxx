#include "WeightCalcCreator.h"
#include "WeightCalcFactory.h"

namespace sbn {
  namespace evwgh {

WeightCalcCreator::WeightCalcCreator(const std::string& wghcalcname) {
  WeightCalcFactory::Register(wghcalcname, this);
}

  }  // namespace evwgh
}  // namespace sbn

