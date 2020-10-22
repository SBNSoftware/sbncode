////////////////////////////////////////////////////////////////////////
// \file    SRTrkSplit.cxx
// \brief   An SRTrkSplit is an object for information from splitting a track
////////////////////////////////////////////////////////////////////////
#include "SRTrkSplit.h"

#include <limits>

namespace caf
{

  SRTrkSplit::SRTrkSplit():
    split(false),
    locAtSplit(std::numeric_limits<double>::signaling_NaN(), std::numeric_limits<double>::signaling_NaN(), std::numeric_limits<double>::signaling_NaN()),
    dirAtSplit(std::numeric_limits<double>::signaling_NaN(), std::numeric_limits<double>::signaling_NaN(), std::numeric_limits<double>::signaling_NaN())
  {  }

  SRTrkSplit::~SRTrkSplit(){  }

  void SRTrkSplit::setDefault()
  {
  }

} // end namespace caf
////////////////////////////////////////////////////////////////////////
