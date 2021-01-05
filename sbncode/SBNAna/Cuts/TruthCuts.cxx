#include "SBNAna/Cuts/TruthCuts.h"
#include "SBNAna/Vars/TruthVars.h"
#include "SBNAna/Cuts/VolumeDefinitions.h"

#include "StandardRecord/Proxy/SRProxy.h"

namespace ana{

  const Cut kTrueActiveVolumeND(
    [](const caf::SRSliceProxy* slc){
      return kHasNu(slc) && PtInVolAbsX(slc->truth.position, avnd);
    }
    );

  const Cut kTrueFiducialVolumeND(
    [](const caf::SRSliceProxy* slc){
      return kHasNu(slc) && PtInVolAbsX(slc->truth.position, fvndAbs);
    }
    );

  const Cut kTrueActiveVolumeFDCryo1(
    [](const caf::SRSliceProxy* slc){
      return kHasNu(slc) && PtInVol(slc->truth.position, avfd_cryo1);
    }
    );

  const Cut kTrueActiveVolumeFDCryo2(
    [](const caf::SRSliceProxy* slc){
      return kHasNu(slc) && PtInVol(slc->truth.position, avfd_cryo2);
    }
    );

  const Cut kIsAntiNu([](const caf::SRSliceProxy* slc){
            if(slc->truth.index < 0) return false;
            return slc->truth.pdg < 0;
  });

  const Cut kIsNu([](const caf::SRSliceProxy* slc){
            if(slc->truth.index < 0) return false;
            return slc->truth.pdg > 0;
  });

  const Cut kHasNu([](const caf::SRSliceProxy* slc){
            return slc->truth.index >= 0;
  });


  const Cut kIsNue([](const caf::SRSliceProxy* slc){
            return slc->truth.index >= 0 && abs(slc->truth.pdg) == 12;
  });

  const Cut kIsNumu([](const caf::SRSliceProxy* slc){
            return slc->truth.index >= 0 && abs(slc->truth.pdg) == 14;
  });

  const Cut kIsNutau([](const caf::SRSliceProxy* slc){
            return slc->truth.index >= 0 && abs(slc->truth.pdg) == 16;
  });


  const Cut kIsNC([](const caf::SRSliceProxy* slc){
            if(slc->truth.index < 0) return false;
            return !slc->truth.iscc;
  });

  const Cut kVtxDistMagCut([](const caf::SRSliceProxy* slc){
            if(slc->truth.index < 0) return true;
            return (kTruthVtxDistMag(slc) < 1);
  });

  const Cut kSlcCompletenessCut([](const caf::SRSliceProxy* slc){
            if(slc->truth.index < 0) return false;
            return (kCompletness(slc) > 0.5);
  });
}
