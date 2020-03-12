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
    density(-5.0),
    len(-5.0),
    open_angle(-5.0),
    dEdx(),
    energy(),
    dir(-5.0, -5.0, -5.0),
    start(-5.0, -5.0, -5.0)
  {
  }

} // end namespace caf
////////////////////////////////////////////////////////////////////////
