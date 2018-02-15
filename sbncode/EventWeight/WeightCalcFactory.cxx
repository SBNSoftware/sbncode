#include "WeightCalcFactory.h" 

namespace sbncode {
namespace evwgh {

  WeightCalc* WeightCalcFactory::Create(const std::string& wghcalcname)
  {
    std::map<std::string, WeightCalcCreator*>::iterator i;
    i = GetTable().find(wghcalcname);
    
    if (i != GetTable().end())
      return i->second->Create();
    else
      return (WeightCalc*)NULL;
  }
  
  void WeightCalcFactory::Register(const std::string& wghcalcname, WeightCalcCreator* creator)
  {
    GetTable()[wghcalcname] = creator;
  }
  
  std::map<std::string, WeightCalcCreator*>& WeightCalcFactory::GetTable()
  {
    static std::map<std::string, WeightCalcCreator*> table;
    return table;
  }

}
}
