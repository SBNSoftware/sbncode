////////////////////////////////////////////////////////////////////////
// \file    SRTruthMatch.cxx
// \brief   SRTruthMatch object for slice-to-neutrino information.
// \author  $Author: psihas@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "SRTruthMatch.h"
#include <limits>

namespace caf
{
  SRTruthMatch::SRTruthMatch():
    visEinslc(std::numeric_limits<float>::signaling_NaN()),
    visEcosmic(std::numeric_limits<float>::signaling_NaN()),
    eff(std::numeric_limits<float>::signaling_NaN()),
    pur(std::numeric_limits<float>::signaling_NaN()),
    index(std::numeric_limits<int>::signaling_NaN()),
    is_numucc_primary(std::numeric_limits<bool>::signaling_NaN())
  {  }

  SRTruthMatch::~SRTruthMatch(){  }

  void SRTruthMatch::setDefault()
  {
    index             = -5;
    visEinslc         = -5.0;
    visEcosmic        = -5.0;
    eff               = -5.0;
    pur               = -5.0;
  }

} // end namespace caf
////////////////////////////////////////////////////////////////////////
