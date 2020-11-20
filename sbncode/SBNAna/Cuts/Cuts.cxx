#include "SBNAna/Cuts/Cuts.h"
#include "SBNAna/Vars/Vars.h"
#include "SBNAna/Cuts/VolumeDefinitions.h"

#include "StandardRecord/Proxy/SRProxy.h"

#include <map>
#include <string>


namespace ana{

  // Dummy cut example for testing
  const SpillCut kFirstEvents = kEvt < 10;

  const SpillCut kFlashTrigger(
    [](const caf::SRSpillProxy* sr){
      return ( sr->pass_flashtrig );
    }
    );

  const SpillCut kCRTHitVetoND(
      [](const caf::SRSpillProxy* sr){
        for (auto const& crtHit: sr->crt_hits){
          if (crtHit.time > -0.1 && crtHit.time < 1.7 && crtHit.position.y>-350 && crtHit.pe>100)
            return false;
        }
        return true;
      }
      );

  const Cut kFiducialVolumeND(
  	[](const caf::SRSliceProxy* slc){
          return PtInVolAbsX(slc->vertex, fvnd);
  	}
  	);

  const Cut kActiveVolumeND(
  	[](const caf::SRSliceProxy* slc){
          return PtInVolAbsX(slc->vertex, avnd);
  	}
  	);

  const Cut kActiveVolumeFDCryo1(
  	[](const caf::SRSliceProxy* slc){
          return PtInVol(slc->vertex, avfd_cryo1);
  	}
  	);

  const Cut kActiveVolumeFDCryo2(
        [](const caf::SRSliceProxy* slc){
          return PtInVol(slc->vertex, avfd_cryo2);
  	}
  	);

} // namespace
