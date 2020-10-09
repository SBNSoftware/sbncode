////////////////////////////////////////////////////////////////////////
// \file    SRTrkRange.cxx
// \brief   An SRTrkRange is an object for p estimates using range
////////////////////////////////////////////////////////////////////////
#include "SRTrkRange.h"

#include <limits>

namespace caf
{

  SRTrkRange::SRTrkRange():
    p_muon(std::numeric_limits<double>::signaling_NaN()),
    p_proton(std::numeric_limits<double>::signaling_NaN())
  {  }

  SRTrkRange::~SRTrkRange(){  }

  void SRTrkRange::setDefault()
  {
    p_muon     = -5.0;
    p_proton   = -5.0;
  }

} // end namespace caf
////////////////////////////////////////////////////////////////////////
