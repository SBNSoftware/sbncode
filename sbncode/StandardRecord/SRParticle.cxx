////////////////////////////////////////////////////////////////////////
// \file    SRParticle.cxx
// \brief   An SRParticle is a high level true track object.  It knows
//          true id, direction, length, but does no hit information.
////////////////////////////////////////////////////////////////////////
#include "SRParticle.h"

namespace caf
{

  SRParticle::SRParticle()
  {
  }

  bool SRParticle::IsPrimary() const {
    return start_process ==  kPrimary;
  }

  bool SRParticle::HasBraggPeak() const {
    return \
      // check contained, id, & end process (stopping-only)
      is_contained &&
      ( abs(pdg) == 13 || abs(pdg) == 2212 || 
	abs(pdg) == 211 || abs(pdg) == 321 ) &&
      ( end_process ==  kCoupledTransportation ||
	end_process ==  kFastScintillation ||
	end_process ==  kDecay ||
	end_process ==  kmuMinusCaptureAtRest );
  }

  bool SRParticle::IsGenie() const {
    return gstatus !=  kNotGenie;
  }

  bool SRParticle::IsStable() const {
    return \
       gstatus ==  kNotGenie // non-genie particles are stable
    || gstatus ==  kIStStableFinalState; // stable genie particle

  }

  

} // end namespace caf
////////////////////////////////////////////////////////////////////////
