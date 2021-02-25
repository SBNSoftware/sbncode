#ifndef _SBN_WEIGHTCALCFACTORY_H_
#define _SBN_WEIGHTCALCFACTORY_H_

// Imported from LArSoft's larsim EventWeight

#include <map>
#include <string>

namespace sbn {
  namespace evwgh {

class WeightCalc;
class WeightCalcCreator;

class WeightCalcFactory {
public:
  static WeightCalc* Create(const std::string& classname);

  static void Register(const std::string& wghcalcname,
                       WeightCalcCreator* creator);

private:
  static std::map<std::string, WeightCalcCreator*>& GetTable();
};

  }  // namespace evwgh
}  // namespace sbn

#endif // _SBN_WEIGHTCALCFACTORY_H_

