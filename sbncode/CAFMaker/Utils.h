//////////////////////////////////////////////////////////////////////
// \file    Utils.h
// \brief   CAF Utils ported from NOvA CAFMaker 
// \author  $Author: psihas@fnal.gov
//////////////////////////////////////////////////////////////////////

#ifndef CAF_UTILS_H
#define CAF_UTILS_H

#include "canvas/Persistency/Common/PtrVector.h"

#include <vector>

namespace caf
{
  template<class T> std::vector<T> PtrVecToVec(const art::PtrVector<T>& xs)
  {
    std::vector<T> ret;
    ret.reserve(xs.size());
    for(const art::Ptr<T>& x: xs) ret.push_back(*x);
    return ret;
  }
}

#endif
