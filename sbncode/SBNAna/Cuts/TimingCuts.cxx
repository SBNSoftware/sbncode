#include "SBNAna/Cuts/TimingCuts.h"

#include "StandardRecord/Proxy/SRProxy.h"

/// Beam spill duration is 1.6 microseconds
/// TO DO: Change numbers to reflect reality:
/// Where the beam spill actually starts after the beam trigger?

//currently unused, 23 Nov 2020
//const double kBeamWindowMinMicroSec = 0;
//const double kBeamWindowMaxMicroSec = 1.6;

namespace ana
{

  // const SpillCut kInBeamSpill(
  // 	[](const caf::SRSliceProxy* slc)
  //   {
  //   	const double t = slc->meantime;
  //   	return (t >= kBeamWindowMinMicroSec && t <= kBeamWindowMaxMicroSec);
  //   }
  //   );

}
