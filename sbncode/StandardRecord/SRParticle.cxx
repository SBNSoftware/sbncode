////////////////////////////////////////////////////////////////////////
// \file    SRParticle.cxx
// \brief   An SRParticle is a high level track object.  It knows its
//          direction and length, but does not own its cell hits.
////////////////////////////////////////////////////////////////////////
#include "SRParticle.h"

namespace caf
{

  SRParticle::SRParticle()
  {
  }

  bool SRParticle::IsPrimary() const {
    return start_process == SRParticle::primary;
  }

  bool SRParticle::HasBraggPeak() const {
    return \
      // check contained
      is_contained &&
      // check particle ID
      (abs(pdg) == 13 || abs(pdg) == 2212 || abs(pdg) == 211 || abs(pdg) == 321) &&
      // check end process (stopping-only)
      (end_process == SRParticle::CoupledTransportation ||
       end_process == SRParticle::FastScintillation ||
       end_process == SRParticle::Decay ||
       end_process == SRParticle::muMinusCaptureAtRest);
  }

  bool SRParticle::IsGenie() const {
    return gstatus != SRParticle::kNotGenie;
  }

  bool SRParticle::IsStable() const {
    return \
       gstatus == SRParticle::kNotGenie // non-genie particles are stable
    || gstatus == SRParticle::kIStStableFinalState; // stable genie particle

  }

  

} // end namespace caf
////////////////////////////////////////////////////////////////////////
