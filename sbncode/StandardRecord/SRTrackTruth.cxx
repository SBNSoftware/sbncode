////////////////////////////////////////////////////////////////////////
// \file    SRTrack.cxx
// \brief   An SRTrack is a high level track object.  It knows its
//          direction and length, but does not own its cell hits.
////////////////////////////////////////////////////////////////////////
#include "SRTrackTruth.h"
#include <bits/stdc++.h>

namespace caf
{

SRTrackTruth::SRTrackTruth():
  total_deposited_energy(std::numeric_limits<float>::signaling_NaN()),
  nmatches(0)
{}

ParticleMatch::ParticleMatch():
  G4ID(INT_MIN),
  energy(std::numeric_limits<float>::signaling_NaN())
{}

} // end namespace caf
////////////////////////////////////////////////////////////////////////
