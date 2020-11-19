////////////////////////////////////////////////////////////////////////
// \file    SRFakeReco.cxx
// \brief   SRFakeReco holds true interaction info.
// \author  $Author: grayputnam@uchicago.edu
////////////////////////////////////////////////////////////////////////

#include "sbncode/StandardRecord/SRFakeReco.h"


namespace caf
{

  SRFakeReco::SRFakeReco():
    nuE(std::numeric_limits<float>::signaling_NaN()),
    wgt(std::numeric_limits<float>::signaling_NaN())
  {  }

} // end namespace caf
////////////////////////////////////////////////////////////////////////
