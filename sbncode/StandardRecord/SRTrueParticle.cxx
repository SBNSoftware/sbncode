////////////////////////////////////////////////////////////////////////
// \file    SRTrueParticle.cxx
// \brief   An SRTrueParticle is a high level true track object. It
//          knows true id, direction, length, but  no hit information.
////////////////////////////////////////////////////////////////////////
#include "sbncode/StandardRecord/SRTrueParticle.h"
#include <climits>

namespace caf
{

  SRTrueParticle::SRTrueParticle():
    planeVisE(std::numeric_limits<float>::signaling_NaN()),
    genE(std::numeric_limits<float>::signaling_NaN()),
    startE(std::numeric_limits<float>::signaling_NaN()),
    endE(std::numeric_limits<float>::signaling_NaN()),
    genT(std::numeric_limits<float>::signaling_NaN()),
    startT(std::numeric_limits<float>::signaling_NaN()),
    endT(std::numeric_limits<float>::signaling_NaN()),
    length(std::numeric_limits<float>::signaling_NaN()),
    wallin(kWallNone),
    wallout(kWallNone),
    cont_tpc(false),
    crosses_tpc(false),
    contained(false),
    pdg(INT_MIN),
    G4ID(INT_MIN),
    interaction_id(INT_MIN),
    generator(kUnknownGenerator),
    start_process(g4_process_(-1)), // TODO do we need an "unknown" process?
    end_process(g4_process_(-1)),
    gstatus(kIStUndefined)
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
