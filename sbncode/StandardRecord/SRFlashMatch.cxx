////////////////////////////////////////////////////////////////////////
// \file    SRFlashMatch.cxx
// \brief   SRFlashMatch object for flashmatch summary information.
// \author  $Author: psihas@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "SRFlashMatch.h"
#include <limits>

namespace caf
{
  SRFlashMatch::SRFlashMatch():
    score(std::numeric_limits<float>::signaling_NaN()),
    time(std::numeric_limits<float>::signaling_NaN()),
    pe(std::numeric_limits<float>::signaling_NaN())
  {  }


  SRFlashMatch::~SRFlashMatch(){  }


  void SRFlashMatch::setDefault()
  {
    present        = false;
    score          = -5.0;
    time           = -5.0;
    pe             = -5.0;
  }

} // end namespace caf
////////////////////////////////////////////////////////////////////////
