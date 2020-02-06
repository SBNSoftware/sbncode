////////////////////////////////////////////////////////////////////////
// \file    SRTrack.cxx
// \brief   An SRTrack is a high level track object.  It knows its
//          direction and length, but does not own its cell hits.
////////////////////////////////////////////////////////////////////////
#include "SRTrackTruth.h"

namespace caf
{
  int SRTrackTruth::GetPrimaryMatchID() const {
    if (matches.size() == 0) return -1;
    return matches[0].G4ID;
  }

  float SRTrackTruth::Purity() const {
    if (matches.size() == 0) return 0.;
    if (total_deposited_energy < 1e-6) return 1.;
    return matches[0].energy / total_deposited_energy;
  }

} // end namespace caf
////////////////////////////////////////////////////////////////////////
