////////////////////////////////////////////////////////////////////////
// \file    SRShower.cxx
// \brief   An SRShower is a high level shower object.  It knows its
//          direction and length, but does not own its cell hits.
////////////////////////////////////////////////////////////////////////
#include "SRShower.h"

namespace caf
{

  SRShower::SRShower():
    bestplane(-5),
    dEdx(-5.0),
    energy(-5.0),
    len(-5.0),
    openAngle(-5.0),
    dir(-5.0, -5.0, -5.0),
    start(-5.0, -5.0, -5.0)
  {
  }

} // end namespace caf
////////////////////////////////////////////////////////////////////////
