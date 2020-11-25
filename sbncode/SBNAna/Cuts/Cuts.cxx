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
          if (crtHit.time > 0. && crtHit.time < 1.8 && crtHit.position.y>-350 && crtHit.pe>100)
            return false;
        }
        return true;
      }
      );

  const SpillCut kCRTHitVetoFD(
      [](const caf::SRSpillProxy* sr){
        for (auto const& crtHit: sr->crt_hits){
          if (crtHit.time > 0. && crtHit.time < 1.8 && crtHit.pe>100) // no discriminating in position, yet.
            return false;
        }
        return true;
      }
      );

  const Cut kFiducialVolumeND(
    [](const caf::SRSliceProxy* slc){
          return PtInVolAbsX(slc->vertex, fvndAbs);
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

  const Cut kSlcIsRecoNu([](const caf::SRSliceProxy *slc)
       {
         return !slc->is_clear_cosmic;
       });

  const Cut kSlcNuScoreCut([](const caf::SRSliceProxy *slc)
       {
         return (kSlcIsRecoNu(slc) && slc->nu_score>0.4);
       });

  const Cut kSlcHasFlashMatch([](const caf::SRSliceProxy *slc)
       {
         return slc->fmatch.present;
       });
  const Cut kSlcFlashMatchCut([](const caf::SRSliceProxy *slc)
       {
         return (kSlcHasFlashMatch(slc) && slc->fmatch.score>0 && slc->fmatch.score<6);
       });
} // namespace
