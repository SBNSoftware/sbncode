////////////////////////////////////////////////////////////////////////
// \file    SRFakeRecoParticle.cxx
// \brief   SRFakeRecoParticle holds true interaction info.
// \author  $Author: grayputnam@uchicago.edu
////////////////////////////////////////////////////////////////////////

#include "SRFakeRecoParticle.h"


namespace caf
{

  SRFakeRecoParticle::SRFakeRecoParticle():
    ke(std::numeric_limits<float>::signaling_NaN()),
    costh(std::numeric_limits<float>::signaling_NaN()),
    len(std::numeric_limits<float>::signaling_NaN()),
    pid(std::numeric_limits<int>::signaling_NaN())
  {  }

} // end namespace caf
////////////////////////////////////////////////////////////////////////
