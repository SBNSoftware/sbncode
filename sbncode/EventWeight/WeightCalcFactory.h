#ifndef _WEIGHTCALCFACTORY_H_
#define _WEIGHTCALCFACTORY_H_

#include "WeightCalc.h"
#include "WeightCalcCreator.h"

namespace sbncode {
namespace evwgh {
  class WeightCalcFactory
  {
  public:
    static WeightCalc* Create(const std::string& classname);
    static void Register(const std::string& wghcalcname, 
			 WeightCalcCreator* creator);
  
  private:
    static std::map<std::string, WeightCalcCreator*>& GetTable();
  };
}
}

#endif // _WEIGHTCALCFACTORY_H_
