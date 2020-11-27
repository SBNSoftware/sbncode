////////////////////////////////////////////////////////////////////////
// \file    SRSlice.cxx
// \brief   SRSlice object for slice summary information.
// \author  $Author: psihas@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "sbncode/StandardRecord/SRSlice.h"

#include <limits>
#include <climits>

namespace caf
{
  FlashMatch::FlashMatch():
    present(false),
    score(std::numeric_limits<float>::signaling_NaN()),
    time(std::numeric_limits<float>::signaling_NaN()),
    pe(std::numeric_limits<float>::signaling_NaN())
  {}

  SRSlice::SRSlice():
    producer(UINT_MAX),
    charge(std::numeric_limits<float>::signaling_NaN()),
    is_clear_cosmic(false),
    nu_pdg(INT_MIN),
    nu_score(std::numeric_limits<float>::signaling_NaN()),
    self(INT_MIN)
  {  }


  SRSlice::~SRSlice(){  }


  void SRSlice::setDefault()
  {
    charge         = -5;
  }


} // end namespace caf
////////////////////////////////////////////////////////////////////////
