#ifndef _sbnanalysis_XGBoostPID_hh_
#define _sbnanalysis_XGBoostPID_hh_

#include <map>
#include <string>

#include "../xgboost/include/xgboost/c_api.h"

namespace ana {
namespace SBNOsc {

/**
 * Class to handle particle ID algorithm using XGBoost BDT
 */
class XGBoostPID {

public:
  XGBoostPID();
  float PredictOne(const std::map<std::string, float> &data) const;
  void SetModelFile(const char *modelFile);
  bool Ready() const { return fReady; }
  ~XGBoostPID();

private:
  BoosterHandle fHandle;
  bool fReady;

};
} // end namespace ana
} // end namespace SBNOsc
#endif
