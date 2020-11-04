#ifndef BASEFLASHHYPOTHESIS_CXX
#define BASEFLASHHYPOTHESIS_CXX

#include "BaseFlashHypothesis.h"

namespace flashmatch {

  Flash_t BaseFlashHypothesis::GetEstimate(const QCluster_t& tpc) const
  {
    Flash_t res;
    //res.pe_v.resize(OpDetXArray().size());

    FillEstimate(tpc,res);
    return res;
  }

}
#endif
