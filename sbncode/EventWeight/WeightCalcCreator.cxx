#include "WeightCalcCreator.h" 
#include "WeightCalcFactory.h"

namespace sbncode {
namespace evwgh {
  WeightCalcCreator::WeightCalcCreator(const std::string& wghcalcname) 
  {
    WeightCalcFactory::Register(wghcalcname, this);
  }
}
}
