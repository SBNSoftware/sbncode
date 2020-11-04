#ifndef OPT0FINDER_BASEFLASHMATCH_CXX
#define OPT0FINDER_BASEFLASHMATCH_CXX

#include "BaseFlashMatch.h"

namespace flashmatch {
  
  Flash_t BaseFlashMatch::GetEstimate(const QCluster_t& tpc) const
  {
    return _flash_hypothesis->GetEstimate(tpc);
  }

  void BaseFlashMatch::FillEstimate(const QCluster_t& tpc, Flash_t& opdet) const
  {
    _flash_hypothesis->FillEstimate(tpc,opdet);
  }

  void BaseFlashMatch::SetFlashHypothesis(flashmatch::BaseFlashHypothesis* alg)
  {
    _flash_hypothesis = alg;
  }

}

#endif
