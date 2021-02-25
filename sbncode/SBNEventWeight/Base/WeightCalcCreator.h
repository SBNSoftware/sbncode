#ifndef _SBN_WEIGHTCALCCREATOR_H_
#define _SBN_WEIGHTCALCCREATOR_H_

// Imported from LArSoft's larsim EventWeight

#include <iostream>
#include <string>

namespace sbn {
  namespace evwgh {

class WeightCalc;

class WeightCalcCreator {
public:
  WeightCalcCreator(const std::string& classname);

  virtual ~WeightCalcCreator() = default;

  virtual WeightCalc* Create() = 0;
};


template <class T>
class WeightCalcImpl : public WeightCalcCreator {
public:
  WeightCalcImpl<T>(const std::string& classname)
      : WeightCalcCreator(classname) {}

  virtual ~WeightCalcImpl<T>() {}

  virtual WeightCalc* Create() { return new T; }
};

  }  // namespace evwgh
}  // namespace sbn


#define DECLARE_WEIGHTCALC(wghcalc)		\
private:					\
  static const sbn::evwgh::WeightCalcImpl<wghcalc> creator;

#define REGISTER_WEIGHTCALC(wghcalc)  		\
const sbn::evwgh::WeightCalcImpl<wghcalc> wghcalc::creator(#wghcalc);

#endif // _SBN_WEIGHTCALCFACTORY_H_

