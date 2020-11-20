#pragma once

// Definition of the generic Cut object
#include "CAFAna/Core/Cut.h"

namespace ana{

  // Select beam mode
  // extern const SpillCut kIsRHC;
  extern const SpillCut kFirstEvents;
  extern const SpillCut kFlashTrigger;

  extern const SpillCut kCRTHitVetoND;

  extern const Cut kActiveVolumeND;
  extern const Cut kActiveVolumeFDCryo1;
  extern const Cut kActiveVolumeFDCryo2;

} // namespace
