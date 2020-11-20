#include "SBNAna/Cuts/TruthCuts.h"
#include "SBNAna/Cuts/VolumeDefinitions.h"

#include "StandardRecord/Proxy/SRProxy.h"

namespace ana{

  const Cut kTrueActiveVolumeND(
    [](const caf::SRSliceProxy* slc){
      if (!kHasNu(slc)) return false;
      bool x = (avnd.xmin < abs(slc->truth.position.x)) && (abs(slc->truth.position.x) < avnd.xmax);
      bool y = (avnd.ymin < slc->truth.position.y) && (slc->truth.position.y < avnd.ymax);
      bool z = (avnd.zmin < slc->truth.position.z) && (slc->truth.position.z < avnd.zmax);
      return(x && y && z);
    }
    );

  const Cut kTrueActiveVolumeFDCryo1(
    [](const caf::SRSliceProxy* slc){
      if (!kHasNu(slc)) return false;
      bool x = (avfd_cryo1.xmin < slc->truth.position.x) && (slc->truth.position.x < avfd_cryo1.xmax);
      bool y = (avfd_cryo1.ymin < slc->truth.position.y) && (slc->truth.position.y < avfd_cryo1.ymax);
      bool z = (avfd_cryo1.zmin < slc->truth.position.z) && (slc->truth.position.z < avfd_cryo1.zmax);
      return(x && y && z);
    }
    );

  const Cut kTrueActiveVolumeFDCryo2(
    [](const caf::SRSliceProxy* slc){
      if (!kHasNu(slc)) return false;
      bool x = (avfd_cryo2.xmin < slc->truth.position.x) && (slc->truth.position.x < avfd_cryo2.xmax);
      bool y = (avfd_cryo2.ymin < slc->truth.position.y) && (slc->truth.position.y < avfd_cryo2.ymax);
      bool z = (avfd_cryo2.zmin < slc->truth.position.z) && (slc->truth.position.z < avfd_cryo2.zmax);
      return(x && y && z);
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


}
