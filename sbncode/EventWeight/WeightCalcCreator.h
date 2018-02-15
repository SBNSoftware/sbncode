#ifndef _WEIGHTCALCCREATOR_H_
#define _WEIGHTCALCCREATOR_H_

#include "WeightCalc.h"

namespace sbncode {
namespace evwgh {
  class WeightCalcCreator
  {
  public:
    WeightCalcCreator(const std::string& classname);
    virtual ~WeightCalcCreator() {};

    virtual WeightCalc* Create() = 0;
  };
  
  template <class T>
    class WeightCalcImpl : public WeightCalcCreator
  {
  public:
    WeightCalcImpl<T>(const std::string& classname) : WeightCalcCreator(classname) {}
    virtual ~WeightCalcImpl<T>() {}
    virtual WeightCalc* Create() { return new T; }
  };

#define DECLARE_WEIGHTCALC(wghcalc)		\
  private:					\
    static const WeightCalcImpl<wghcalc> creator;

#define REGISTER_WEIGHTCALC(wghcalc)  		\
  const WeightCalcImpl<wghcalc> wghcalc::creator(#wghcalc);
}
}

#endif // _WEIGHTCALCFACTORY_H_
