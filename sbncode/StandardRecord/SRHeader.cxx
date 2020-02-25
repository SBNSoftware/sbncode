////////////////////////////////////////////////////////////////////////
// \file    SRHeader.cxx
// \brief   SRHeader object for slice id information.
// \author  $Author: psihas@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "SRHeader.h"

#include <limits>


namespace caf
{
  SRHeader::SRHeader():
  run(0),
  subrun(0),
  // cycle(0),
  // batch(0),
  evt(0),
  subevt(0),
  ismc(false),
  det(std::numeric_limits<Det_t>::signaling_NaN())
  // blind(false),
  {  }

  SRHeader::~SRHeader(){  }

} // end namespace caf
////////////////////////////////////////////////////////////////////////
