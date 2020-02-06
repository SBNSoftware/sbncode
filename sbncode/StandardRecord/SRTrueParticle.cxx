////////////////////////////////////////////////////////////////////////
// \file    SRTrueParticle.cxx
// \brief   An SRTrueParticle is a high level true track object. It
//          knows true id, direction, length, but  no hit information.
////////////////////////////////////////////////////////////////////////
#include "SRTrueParticle.h"

namespace caf
{

  SRTrueParticle::SRTrueParticle()
  {
  }

  bool SRTrueParticle::IsPrimary() const {
    return start_process == kG4primary;
  }

  bool SRTrueParticle::HasBraggPeak() const {
    return \
      // check contained, id, & end process (stopping-only)
      contained &&
      ( abs(pdg) == 13 || abs(pdg) == 2212 || 
	abs(pdg) == 211 || abs(pdg) == 321 ) &&
      ( end_process ==  kG4CoupledTransportation ||
	end_process ==  kG4FastScintillation ||
	end_process ==  kG4Decay ||
	end_process ==  kG4muMinusCaptureAtRest );
  }

  bool SRTrueParticle::IsGenie() const {
    return gstatus !=  kNotGenie;
  }

  bool SRTrueParticle::IsStable() const {
    return \
       gstatus ==  kNotGenie // non-genie particles are stable
    || gstatus ==  kIStStableFinalState; // stable genie particle

  }

} // end namespace caf
////////////////////////////////////////////////////////////////////////
