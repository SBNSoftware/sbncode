////////////////////////////////////////////////////////////////////////
// \file    SRTrackCalo.cxx
// \brief   An SRTrackCalo is a high level track ParticlePID object for
//          for larana Chi2ParticleID results. 
////////////////////////////////////////////////////////////////////////

#include "sbncode/StandardRecord/SRTrackCalo.h"

#include <limits>

namespace caf
{

  SRTrackCalo::SRTrackCalo():
    nhit(std::numeric_limits<int>::signaling_NaN()),
    ke(std::numeric_limits<double>::signaling_NaN()),
    charge(std::numeric_limits<double>::signaling_NaN())
  {  }

  SRTrackCalo::~SRTrackCalo(){  }

  void SRTrackCalo::setDefault()
  {
    nhit          = -1;
    ke            = -1.;
    charge        = -1.;
  }

} // end namespace caf
////////////////////////////////////////////////////////////////////////
