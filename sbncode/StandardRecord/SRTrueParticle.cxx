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
    return start_process ==  kPrimary;
  }

  bool SRTrueParticle::HasBraggPeak() const {
    return \
      // check contained, id, & end process (stopping-only)
      isContained &&
      ( abs(pdg) == 13 || abs(pdg) == 2212 || 
	abs(pdg) == 211 || abs(pdg) == 321 ) &&
      ( end_process ==  kCoupledTransportation ||
	end_process ==  kFastScintillation ||
	end_process ==  kDecay ||
	end_process ==  kMuMinusCaptureAtRest );
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
