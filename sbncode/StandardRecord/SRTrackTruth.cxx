////////////////////////////////////////////////////////////////////////
// \file    SRTrack.cxx
// \brief   An SRTrack is a high level track object.  It knows its
//          direction and length, but does not own its cell hits.
////////////////////////////////////////////////////////////////////////
#include "SRTrackTruth.h"

namespace caf
{

SRTrackTruth::SRTrackTruth():
  total_deposited_energy(-1.),
  nmatches(0)
{}

ParticleMatch::ParticleMatch():
  G4ID(-1.),
  energy(-1.)
{}

} // end namespace caf
////////////////////////////////////////////////////////////////////////
