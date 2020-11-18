////////////////////////////////////////////////////////////////////////
// \file    SRTrkRange.cxx
// \brief   An SRTrkRange is an object for p estimates using range
////////////////////////////////////////////////////////////////////////

#include "sbncode/StandardRecord/SRTrkRange.h"

#include <limits>

namespace caf
{

  SRTrkRange::SRTrkRange():
    p_muon(std::numeric_limits<float>::signaling_NaN()),
    p_proton(std::numeric_limits<float>::signaling_NaN())
  {  }

  SRTrkRange::~SRTrkRange(){  }

  void SRTrkRange::setDefault()
  {
    p_muon     = -5.0;
    p_proton   = -5.0;
  }

} // end namespace caf
////////////////////////////////////////////////////////////////////////
