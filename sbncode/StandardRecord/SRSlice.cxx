////////////////////////////////////////////////////////////////////////
// \file    SRSlice.cxx
// \brief   SRSlice object for slice summary information.
// \author  $Author: psihas@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "SRSlice.h"
#include <limits>

namespace caf
{
  SRSlice::SRSlice():
    charge(std::numeric_limits<float>::signaling_NaN())
  {  }


  SRSlice::~SRSlice(){  }


  void SRSlice::setDefault()
  {
    charge         = -5;
  }


} // end namespace caf
////////////////////////////////////////////////////////////////////////
